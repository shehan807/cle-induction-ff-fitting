"""Reproduce MBIS partial atomic charges (HF/aug-cc-pVDZ) for four monomers:
acetic acid, formamide, imidazole, pyrazole. Geometries in geometries/.

Env: environment.yml.  Run: python partial_charges.py
"""
import csv
import warnings
from pathlib import Path

import numpy as np
from openff.toolkit import Molecule
from pympfit import (
    MBISSettings,
    MoleculeMBISRecord,
    Psi4MBISGenerator,
    extract_mbis_charges,
)

warnings.filterwarnings("ignore", message=r".*from_pdb.*")

SETTINGS = MBISSettings(
    method="scf", basis="aug-cc-pvdz",
    max_moment=3, max_radial_moment=4, limit=3,
    multipole_format="spherical", mbis_d_convergence=9,
    mbis_radial_points=99, mbis_spherical_points=590,
)

# smiles set bonds/formal charges; geometry comes from the PDB.
MOLECULES = {
    "acetic_acid": "CC(=O)O",
    "formamide": "NC=O",
    "imidazole": "c1cnc[nH]1",
    "pyrazole": "c1cc[nH]n1",
}


def mbis_charges(name, smiles):
    """Psi4 SCF + MBIS on geometries/<name>.pdb -> (elements, labels, charges)."""
    pdb = Path(__file__).resolve().parent / "geometries" / f"{name}.pdb"
    mol = Molecule.from_pdb_and_smiles(str(pdb), smiles)
    conf, mp = Psi4MBISGenerator.generate(mol, mol.conformers[0], SETTINGS)
    rec = MoleculeMBISRecord.from_molecule(mol, conf, mp, SETTINGS)
    q = np.asarray(extract_mbis_charges(rec).value, dtype=float)
    lines = pdb.read_text().splitlines()
    atoms = [ln.split() for ln in lines if ln.startswith(("ATOM", "HETATM"))]
    elements = [a.symbol for a in mol.atoms]
    assert elements == [a[-1] for a in atoms], f"{name}: PDB/openff atom order mismatch"
    return elements, [a[2] for a in atoms], q


def write2csv(rows, path="mbis_charges.csv"):
    """Write (molecule, atom_idx, element, label, charge) rows to `path`."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(("molecule", "atom_idx", "element", "label", "charge"))
        w.writerows(rows)


def test(rows):
    """Every monomer is neutral, and imidazole N00 reproduces the SI value (-0.5362)."""
    chg = {(name, lbl): c for name, _, _, lbl, c in rows}
    for name in MOLECULES:
        net = sum(c for (m, _), c in chg.items() if m == name)
        assert abs(net) < 1e-3, f"{name} not neutral: net {net:+.2e}"
    n00 = chg[("imidazole", "N00")]
    assert f"{n00:+.4f}" == "-0.5362", f"imidazole N00 = {n00:+.6f}"


def main():
    rows = []
    for name, smiles in MOLECULES.items():
        elements, labels, q = mbis_charges(name, smiles)
        print(f"\n{name}  (net {q.sum():+.4f})")
        for i, (e, lbl, c) in enumerate(zip(elements, labels, q)):
            print(f"  {i:2d}  {e:<2s} {lbl:<4s} {c:+.6f}")
            rows.append((name, i, e, lbl, c))
    write2csv(rows)
    test(rows)


if __name__ == "__main__":
    main()
