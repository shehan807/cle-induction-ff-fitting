# MBIS charges

MBIS partial atomic charges (HF/aug-cc-pVDZ) for acetic acid, formamide, imidazole,
and pyrazole, computed on the geometries in `geometries/`.

Setup:

    micromamba env create -f environment.yml

Run:

    micromamba run -n cle-mbis python partial_charges.py

Prints the charges, writes `mbis_charges.csv`, and checks imidazole's N00 against the
reference value. Takes about a minute.

Charges are computed on the PDB coordinates directly; the SMILES only set bond orders
and formal charges, so don't swap in a SMILES-generated geometry.
