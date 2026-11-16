"""Checks on the force field XMLs in force_fields/ (run partial_charges.py, polarizabilities.py and
force_field_xml.py first).  Run: python -m pytest test_force_field_xml.py"""
import csv
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u
import pytest

from force_field_xml import BOHR3_TO_NM3, K, MOLECULES
from partial_charges import METHODS
from polarizabilities import THOLE

HERE = Path(__file__).resolve().parent
TOL = 2e-8                                   # values are written with 8 decimals


def table(path, column):
    rows = {}
    with open(HERE / path) as f:
        for r in csv.DictReader(f):
            rows.setdefault(r["molecule"], {})[r["label"]] = float(r[column])
    return rows


CHARGES = {m.lower(): table("partial_charges.csv", m.lower()) for m in METHODS}
ALPHA_A3 = table("polarizabilities.csv", "alpha_A3")
CASES = [(m, n) for m in CHARGES for n in MOLECULES]          # (charge method, molecule)


def parse(method, name):
    """(charge by type, {atom label: (drude charge, polarizability nm^3, thole)}, label by type)."""
    root = ET.parse(HERE / "force_fields" / method / f"{name}.xml").getroot()
    label_of = {a.get("type"): a.get("name") for a in root.find("Residues/Residue").findall("Atom")}
    charge = {a.get("type"): float(a.get("charge")) for a in root.find("NonbondedForce")}
    drude = {label_of[p.get("type2")]: (float(p.get("charge")), float(p.get("polarizability")), float(p.get("thole")))
             for p in root.find("DrudeForce")}
    return charge, drude, label_of


@pytest.mark.parametrize("method, name", CASES)
def test_net_charge_is_zero(method, name):
    assert abs(sum(parse(method, name)[0].values())) < 1e-6


@pytest.mark.parametrize("method, name", CASES)
def test_core_plus_drude_is_the_input_charge(method, name):
    charge, drude, label_of = parse(method, name)
    total = {}
    for t, q in charge.items():                      # a Drude's label is its atom's label with a "D" prefix
        lab = label_of[t]
        lab = lab[1:] if lab not in CHARGES[method][name] else lab
        total[lab] = total.get(lab, 0.0) + q
    assert total.keys() == CHARGES[method][name].keys()
    assert all(abs(total[lab] - CHARGES[method][name][lab]) < TOL for lab in total)


@pytest.mark.parametrize("method, name", CASES)
def test_drude_charge_and_polarizability_agree(method, name):
    """alpha = q_D^2 / k, and the same q_D sits in NonbondedForce and DrudeForce."""
    charge, drude, label_of = parse(method, name)
    for lab, (q, alpha_nm3, _) in drude.items():
        assert q < 0
        assert alpha_nm3 == pytest.approx(q * q / K * BOHR3_TO_NM3, rel=1e-5)
        assert charge[next(t for t, l in label_of.items() if l == "D" + lab)] == q


@pytest.mark.parametrize("method, name", CASES)
def test_polarizabilities_are_the_fitted_ones_in_sensible_units(method, name):
    """Every heavy atom (and no H) has a Drude; alpha in nm^3 x 1000 = the fitted A^3, 0.5-3 A^3."""
    _, drude, _ = parse(method, name)
    heavy = {lab for lab, a in ALPHA_A3[name].items() if a > 0}
    assert drude.keys() == heavy
    for lab, (_, alpha_nm3, thole) in drude.items():
        assert 1000 * alpha_nm3 == pytest.approx(ALPHA_A3[name][lab], rel=1e-5)
        assert 0.5 < 1000 * alpha_nm3 < 3.0
        assert 2 * thole == THOLE


@pytest.mark.parametrize("method, name", CASES)
def test_only_charges_and_polarizabilities_changed(method, name):
    """Everything else in the XML is the template, unchanged."""
    def strip(path):
        root = ET.parse(path).getroot()
        for el in root.iter():
            el.attrib.pop("polarizability", None)
            if el.tag in ("Atom", "Particle"):
                el.attrib.pop("charge", None)
        return ET.tostring(root)
    assert strip(HERE / "force_fields" / method / f"{name}.xml") == strip(HERE / "templates" / f"{name}.xml")


def openmm_system(method, name, shifts_nm):
    """OpenMM system of one monomer per shift (along x), Drudes added from the XML."""
    ff = app.ForceField(str(HERE / "force_fields" / method / f"{name}.xml"))
    pdb = app.PDBFile(str(HERE / "geometries" / f"{name}.pdb"))
    res = ET.parse(HERE / "templates" / f"{name}.xml").getroot().find("Residues/Residue")
    xyz = np.array(pdb.positions.value_in_unit(u.nanometer))
    top, pos = app.Topology(), []
    for s in shifts_nm:
        r = top.addResidue(res.get("name"), top.addChain())
        atoms = [top.addAtom(a.name, a.element, r) for a in pdb.topology.atoms()]
        for b in res.findall("Bond"):
            top.addBond(atoms[int(b.get("from"))], atoms[int(b.get("to"))])
        pos += [mm.Vec3(*p) for p in xyz + [s, 0, 0]]
    model = app.Modeller(top, pos * u.nanometer)
    model.addExtraParticles(ff)
    return ff.createSystem(model.topology, nonbondedMethod=app.NoCutoff), model.positions


def scf_energy(method, name, shifts_nm):
    """Energy (kJ/mol) with the Drudes relaxed (self-consistent)."""
    system, pos = openmm_system(method, name, shifts_nm)
    integrator = mm.DrudeSCFIntegrator(0.0005 * u.picoseconds)
    integrator.setMinimizationErrorTolerance(1e-8)
    ctx = mm.Context(system, integrator, mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(pos)
    integrator.step(1)
    return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)


@pytest.mark.parametrize("method, name", CASES)
def test_openmm_reads_the_same_numbers(method, name):
    """OpenMM sees a neutral molecule, the XML's polarizabilities (nm^3) and pair Thole = THOLE."""
    system, _ = openmm_system(method, name, [0.0])
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    dr = next(f for f in system.getForces() if isinstance(f, mm.DrudeForce))
    assert abs(sum(nb.getParticleParameters(i)[0].value_in_unit(u.elementary_charge)
                   for i in range(nb.getNumParticles()))) < 1e-6
    loaded = sorted(dr.getParticleParameters(i)[6].value_in_unit(u.nanometer ** 3) for i in range(dr.getNumParticles()))
    assert loaded == pytest.approx(sorted(a for _, a, _ in parse(method, name)[1].values()), rel=1e-6)
    assert all(dr.getScreenedPairParameters(i)[2] == THOLE for i in range(dr.getNumScreenedPairs()))


@pytest.mark.parametrize("method, name", CASES)
def test_far_dimer_is_finite_and_decays_as_dipole_dipole(method, name):
    """Two copies 2 and 4 nm apart: the interaction is finite, small, and falls by ~2^3 = 8."""
    mono = scf_energy(method, name, [0.0])
    e2, e4 = (scf_energy(method, name, [0.0, r]) - 2 * mono for r in (2.0, 4.0))
    assert np.isfinite([e2, e4]).all() and abs(e2) < 1.0
    assert 6 < e2 / e4 < 10


@pytest.mark.parametrize("method, name", CASES)
def test_charges_reproduce_the_qm_dipole(method, name):
    """Independent check on the charge labels: the XML's point charges give the HF/aug-cc-pVDZ
    dipole to within 0.6 D (largest seen 0.58 D; two MBIS heteroatom charges swapped: 0.9-3 D off).
    MPFIT fits charges to local multipoles rather than the ESP, so it gets 1.0 D (pyrazole: 0.94 D)."""
    import psi4
    from drude_fit.inputs import read_geometry
    labels, elements, xyz = read_geometry(HERE / "geometries" / f"{name}.pdb")
    charge, _, label_of = parse(method, name)
    q = {}
    for t, c in charge.items():
        lab = label_of[t]
        q[lab[1:] if lab not in labels else lab] = q.get(lab[1:] if lab not in labels else lab, 0.0) + c
    dipole = sum(q[lab] * r for lab, r in zip(labels, xyz))
    psi4.core.be_quiet()
    mol = psi4.geometry("0 1\n" + "\n".join(f"{e} {x} {y} {z}" for e, (x, y, z) in zip(elements, xyz))
                        + "\nunits bohr\nno_reorient\nno_com\nsymmetry c1")
    _, wfn = psi4.energy("scf/aug-cc-pvdz", return_wfn=True)
    psi4.oeprop(wfn, "DIPOLE")
    limit = 1.0 if method == "mpfit" else 0.6
    assert np.linalg.norm(dipole - np.array(wfn.variable("SCF DIPOLE"))) * 2.541746 < limit   # e bohr -> D
