"""Write OpenMM Drude force field XMLs for acetic acid, formamide, imidazole and pyrazole from the
charges in partial_charges.csv (one force field per charge method; MBIS is the one used in this work)
and the polarizabilities in polarizabilities.csv (run partial_charges.py and polarizabilities.py first).

Each heavy atom carries a Drude particle of charge q_D = -sqrt(k alpha), k = 0.1 au; its core
keeps q - q_D, so the atom's total charge is still q. The templates in templates/
give the atom types, bonds and Thole parameters; only charges and polarizabilities are filled in.

Run: python force_field_xml.py   (writes force_fields/<method>/<molecule>.xml)
"""
import csv
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from drude_fit.inputs import ANG_TO_BOHR
from partial_charges import METHODS
from polarizabilities import THOLE

HERE = Path(__file__).resolve().parent
MOLECULES = ["acetic_acid", "formamide", "imidazole", "pyrazole"]
K = 0.1                                  # Drude spring constant (au)
BOHR3_TO_NM3 = (0.1 / ANG_TO_BOHR) ** 3  # OpenMM wants polarizabilities in nm^3


def read(path, column):
    """{molecule: {atom label: value}} from one column of a CSV."""
    out = {}
    with open(HERE / path) as f:
        for row in csv.DictReader(f):
            out.setdefault(row["molecule"], {})[row["label"]] = float(row[column])
    return out


def write_xml(name, charge, alpha, method="mbis"):
    """Fill templates/<name>.xml with charges (e) and polarizabilities (bohr^3) by atom label."""
    tree = ET.parse(HERE / "templates" / f"{name}.xml")
    root = tree.getroot()
    type_of = {a.get("name"): a.get("type") for a in root.iter("Atom") if a.get("name")}
    label_of = {t: lab for lab, t in type_of.items()}
    nonbonded = {a.get("type"): a for a in root.find("NonbondedForce")}
    for lab, t in type_of.items():                           # charge of every atom (and Drude)
        if t in nonbonded and lab in charge:
            nonbonded[t].set("charge", f"{charge[lab]:.8f}")
    for p in root.find("DrudeForce"):                        # type1 = Drude, type2 = its atom
        # OpenMM screens a Drude pair with the sum of the two thole values: must equal the fit's
        assert 2 * float(p.get("thole")) == THOLE, f"{name}: template thole {p.get('thole')}"
        lab = label_of[p.get("type2")]
        q_drude = -np.sqrt(K * alpha[lab])
        p.set("charge", f"{q_drude:.8f}")
        p.set("polarizability", f"{alpha[lab] * BOHR3_TO_NM3:.8f}")
        nonbonded[p.get("type1")].set("charge", f"{q_drude:.8f}")
        nonbonded[p.get("type2")].set("charge", f"{charge[lab] - q_drude:.8f}")
    out = HERE / "force_fields" / method / f"{name}.xml"
    out.parent.mkdir(parents=True, exist_ok=True)
    tree.write(out)
    return out


def main():
    alphas = read("polarizabilities.csv", "alpha_bohr3")
    for method in (m.lower() for m in METHODS):
        charges = read("partial_charges.csv", method)
        for name in MOLECULES:
            print(f"wrote {write_xml(name, charges[name], alphas[name], method).relative_to(HERE)}")


if __name__ == "__main__":
    main()
