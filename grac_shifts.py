"""Reproduce the GRAC shifts (PBE0/aug-cc-pVTZ) for four monomers:
acetic acid, formamide, imidazole, pyrazole. Geometries in geometries/.

grac = IP + HOMO, with IP = E(cation) - E(neutral) at the fixed PDB geometry.
Psi4 applies the GRAC asymptotic correction (Gruning 2001); its shift is the ionization
potential plus the (negative) HOMO energy (CAMCASP users_guide). CAMCASP uses it for
the PBE0 response the polarizabilities are fitted to.

Env: environment.yml.  Run: python grac_shifts.py
"""
import csv
from pathlib import Path

import psi4

from drude_fit.inputs import read_geometry

HERE = Path(__file__).resolve().parent
MOLECULES = ["acetic_acid", "formamide", "imidazole", "pyrazole"]

# hartree; Psi4 SAPT(DFT) GRAC runs on the same geometries (paper repo, grac_shifts.yaml)
REFERENCE = {
    "acetic_acid": 0.09039351,
    "formamide": 0.09500891,
    "imidazole": 0.08013530,
    "pyrazole": 0.08197448,
}


def grac_shift(elements, xyz_bohr):
    """Return (ip, homo, grac) in atomic units from a neutral + cation PBE0/aug-cc-pVTZ SCF."""
    body = "\n".join(f"{e} {x:.8f} {y:.8f} {z:.8f}" for e, (x, y, z) in zip(elements, xyz_bohr))
    body += "\nunits bohr\nsymmetry c1\n"

    psi4.set_options({"basis": "aug-cc-pvtz", "scf_type": "df", "reference": "rks"})
    e_neutral, wfn = psi4.energy("pbe0", molecule=psi4.geometry("0 1\n" + body), return_wfn=True)
    homo = float(wfn.epsilon_a().to_array()[wfn.nalpha() - 1])

    psi4.set_options({"reference": "uks"})
    e_cation = psi4.energy("pbe0", molecule=psi4.geometry("1 2\n" + body))

    ip = e_cation - e_neutral
    return ip, homo, ip + homo


def write2csv(rows, path="grac_shifts.csv"):
    """Write (molecule, ip, homo, grac) rows to `path`, in hartree."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(("molecule", "ip_au", "homo_au", "grac_shift_au"))
        w.writerows(rows)


def test(rows):
    """Every shift reproduces the reference SAPT(DFT) value to 1e-5 hartree."""
    for name, _, _, grac in rows:
        diff = grac - REFERENCE[name]
        assert abs(diff) < 1e-5, f"{name}: grac {grac:.8f} vs reference {REFERENCE[name]:.8f}"


def main():
    psi4.set_memory("8 GB")
    psi4.set_num_threads(8)
    psi4.core.set_output_file(str(HERE / "grac_shifts.out"), False)
    rows = []
    print(f"{'molecule':12s} {'IP':>11s} {'HOMO':>11s} {'GRAC':>11s} {'reference':>11s} {'diff':>9s}")
    for name in MOLECULES:
        _, elements, xyz = read_geometry(HERE / "geometries" / f"{name}.pdb")
        ip, homo, grac = grac_shift(elements, xyz)
        print(f"{name:12s} {ip:11.8f} {homo:11.8f} {grac:11.8f} {REFERENCE[name]:11.8f} {grac - REFERENCE[name]:9.1e}")
        rows.append((name, ip, homo, grac))
    write2csv(rows)
    test(rows)


if __name__ == "__main__":
    main()
