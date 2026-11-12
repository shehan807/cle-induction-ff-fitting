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

# GRAC shifts

GRAC shift = IP + HOMO (PBE0/aug-cc-pVTZ), from a neutral and a cation SCF on the
same geometries. CAMCASP needs it for the PBE0 response that the polarizabilities
are fitted to.

Run:

    micromamba run -n cle-mbis python grac_shifts.py

Prints IP, HOMO and the shift (hartree), writes `grac_shifts.csv`, and checks each
shift against the reference SAPT(DFT) value to 1e-5. Takes about a minute.

The geometries are the crystal monomers (first dimer, monomer A, of each crystal's
SAPT dimer set), the same ones used for the charges and the CAMCASP runs.
