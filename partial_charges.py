"""Partial atomic charges for acetic acid, formamide, imidazole and pyrazole by four methods:
MBIS (used in this work), RESP, MPFIT (fit to a GDMA multipole expansion) and AM1-BCC.
MBIS, RESP and MPFIT use HF/aug-cc-pVDZ on the geometries in geometries/; AM1-BCC relaxes
the geometry with AM1 first.

Env: environment.yml.  Run: python partial_charges.py
"""
import contextlib
import csv
import tempfile
import warnings
from pathlib import Path

import numpy as np
from openff.toolkit import Molecule

warnings.filterwarnings("ignore", message=r".*from_pdb.*")

BASIS = "aug-cc-pvdz"
METHODS = ["MBIS", "RESP", "MPFIT", "AM1BCC"]

# smiles set bonds/formal charges; geometry comes from the PDB.
MOLECULES = {
    "acetic_acid": "CC(=O)O",
    "formamide": "NC=O",
    "imidazole": "c1cnc[nH]1",
    "pyrazole": "c1cc[nH]n1",
}


def partial_charges(method, mol):
    """Charges (e) of `mol` (openff Molecule with one conformer) in its atom order."""
    xyz = mol.conformers[0]
    if method == "MBIS":        # density partitioned into atoms (pympfit + Psi4)
        from pympfit import MBISSettings, MoleculeMBISRecord, Psi4MBISGenerator, extract_mbis_charges
        settings = MBISSettings(method="scf", basis=BASIS, max_moment=3, max_radial_moment=4, limit=3,
                                multipole_format="spherical", mbis_d_convergence=9,
                                mbis_radial_points=99, mbis_spherical_points=590)
        conf, mp = Psi4MBISGenerator.generate(mol, xyz, settings)
        return extract_mbis_charges(MoleculeMBISRecord.from_molecule(mol, conf, mp, settings)).value

    if method == "MPFIT":       # Stone's analytical DMA to rank 8, then charges fit to those multipoles
        from pympfit import (GDMASettings, MoleculeGDMARecord, MPFITSVDSolver, Psi4GDMAGenerator,
                             generate_mpfit_charge_parameter)
        settings = GDMASettings(method="scf", basis=BASIS, limit=8, switch=0.0,
                                radius=["H", 0.325, "C", 0.65, "N", 0.65, "O", 0.65])
        conf, mp = Psi4GDMAGenerator.generate(mol, xyz, settings)
        record = MoleculeGDMARecord.from_molecule(mol, conf, mp, settings)
        return generate_mpfit_charge_parameter([record], solver=MPFITSVDSolver(svd_threshold=1e-4)).value

    if method == "RESP":        # two-stage restrained fit to the ESP on 1.4-2.0 x vdW shells (Psi4 resp)
        import psi4
        import resp
        atoms = "\n".join(f"{a.symbol} {x:.10f} {y:.10f} {z:.10f}"
                          for a, (x, y, z) in zip(mol.atoms, xyz.m_as("angstrom")))
        p4 = psi4.geometry(f"0 1\n{atoms}\nunits angstrom\nno_reorient\nno_com\nsymmetry c1")
        p4.set_name("mol")
        opts = {"VDW_SCALE_FACTORS": [1.4, 1.6, 1.8, 2.0], "VDW_POINT_DENSITY": 1.0,
                "METHOD_ESP": "scf", "BASIS_ESP": BASIS,
                "RESTRAINT": True, "RESP_A": 0.0005, "RESP_B": 0.1, "IHFREE": True}
        psi4.core.be_quiet()
        with tempfile.TemporaryDirectory() as d, contextlib.chdir(d):     # resp writes to the cwd
            stage1 = resp.resp([p4], opts)
            resp.set_stage2_constraint(p4, stage1[1], opts)
            opts.update(RESP_A=0.001, grid=["1_mol_grid.dat"], esp=["1_mol_grid_esp.dat"])
            return resp.resp([p4], opts)[1]

    if method == "AM1BCC":      # AM1 + bond-charge corrections (AmberTools via openff-recharge)
        from openff.recharge.charges.qc import QCChargeGenerator, QCChargeSettings
        return QCChargeGenerator.generate(mol, [xyz], QCChargeSettings(theory="am1bcc"))

    raise ValueError(method)


def load(name):
    """openff Molecule for geometries/<name>.pdb, and the PDB atom labels."""
    pdb = Path(__file__).resolve().parent / "geometries" / f"{name}.pdb"
    mol = Molecule.from_pdb_and_smiles(str(pdb), MOLECULES[name])
    atoms = [ln.split() for ln in pdb.read_text().splitlines() if ln.startswith(("ATOM", "HETATM"))]
    assert [a.symbol for a in mol.atoms] == [a[-1] for a in atoms], f"{name}: PDB/openff atom order mismatch"
    return mol, [a[2] for a in atoms]


def write2csv(rows, path="partial_charges.csv"):
    """Write (molecule, atom_idx, element, label, one charge per method) rows to `path`."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(("molecule", "atom_idx", "element", "label", *(m.lower() for m in METHODS)))
        w.writerows(rows)


def test(rows):
    """Every monomer is neutral for every method, and imidazole N00 reproduces the SI MBIS value (-0.5362)."""
    for name in MOLECULES:
        for k, m in enumerate(METHODS):
            net = sum(r[4 + k] for r in rows if r[0] == name)
            assert abs(net) < 1e-3, f"{name} {m} not neutral: net {net:+.2e}"
    n00 = next(r[4] for r in rows if r[0] == "imidazole" and r[3] == "N00")
    assert f"{n00:+.4f}" == "-0.5362", f"imidazole N00 = {n00:+.6f}"


def main():
    rows = []
    for name in MOLECULES:
        mol, labels = load(name)
        q = np.column_stack([np.asarray(partial_charges(m, mol), dtype=float).reshape(-1) for m in METHODS])
        print(f"\n{name}" + "".join(f"{m:>10s}" for m in METHODS))
        for i, (a, lab) in enumerate(zip(mol.atoms, labels)):
            print(f"  {i:2d}  {a.symbol:<2s} {lab:<4s}" + "".join(f"{c:+10.4f}" for c in q[i]))
            rows.append((name, i, a.symbol, lab, *q[i]))
    write2csv(rows)
    test(rows)


if __name__ == "__main__":
    main()
