"""Write the Fortran fitter's two inputs for one molecule, one type per heavy atom:
<mol>.clt (geometry in bohr) and <mol>.constraint (starting polarizability per type; 0 = no Drude).
Same format as SAPT-FF-CLE ff_fitting/pol/polfit/io.py:fit_inputs_from_pdb (per_atom=True).

Run from the tutorial root: python legacy_fortran/fortran_inputs.py <molecule> <out_dir>
"""
import sys
from pathlib import Path

sys.path.insert(0, ".")                       # the tutorial root
from drude_fit.inputs import read_geometry
from polarizabilities import START

Z = {"H": 1, "C": 6, "N": 7, "O": 8}

name, out = sys.argv[1], Path(sys.argv[2])
labels, elements, xyz = read_geometry(Path("geometries") / f"{name}.pdb")
atoms = [f"{lab:<4} {Z[el]:>6.1f} {x:>16.8f} {y:>16.8f} {z:>16.8f}   Type {lab}"
         for lab, el, (x, y, z) in zip(labels, elements, xyz)]
(out / f"{name}.clt").write_text("\n".join([f"Molecule {name}", "Charge 0", *atoms, "End", ""]))
lines = ["! Per atom type: a rank-1 line, then the starting polarizability (bohr^3).",
         "! positive = fit, negative = fix, zero = no Drude.", "!"]
for lab, el in zip(labels, elements):
    lines += [lab, "1", str(START.get(el, 0.0))]
(out / f"{name}.constraint").write_text("\n".join(lines) + "\n")
