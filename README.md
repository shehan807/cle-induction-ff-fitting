# Partial charges

Partial atomic charges for acetic acid, formamide, imidazole, and pyrazole by four
methods: MBIS (used in this work), RESP, MPFIT, and AM1-BCC. MBIS, RESP and MPFIT use
HF/aug-cc-pVDZ on the geometries in `geometries/`; AM1-BCC relaxes the geometry first.

Setup:

    micromamba env create -f environment.yml

Run:

    micromamba run -n cle-mbis python partial_charges.py

Prints one column of charges per method, writes `partial_charges.csv`, and checks
imidazole's MBIS N00 against the reference value.

Charges are computed on the PDB coordinates directly; the SMILES only set bond orders
and formal charges, so don't swap in a SMILES-generated geometry.

# GRAC shifts

GRAC shift = IP + HOMO (PBE0/aug-cc-pVTZ), from a neutral and a cation SCF on the
same geometries. CamCASP needs it for the PBE0 response that the polarizabilities
are fitted to.

Run:

    micromamba run -n cle-mbis python grac_shifts.py

Prints IP, HOMO and the shift (hartree), writes `grac_shifts.csv`, and checks each
shift against the reference SAPT(DFT) value to 1e-5. Takes about a minute.

The geometries are the crystal monomers (first dimer, monomer A, of each crystal's
SAPT dimer set), the same ones used for the charges and the CamCASP runs.

# CamCASP responses

The point-to-point (p2p) response of each monomer (PBE0/d-aug-cc-pVTZ, ALDA+CHF kernel,
GRAC shift from `grac_shifts.csv`) on 2000 random grid points. The polarizabilities are
fitted to these.

Needs CamCASP 6.1 (https://gitlab.com/anthonyjs/camcasp) with `$CAMCASP` set to its
install; keep that path under 80 characters, since CamCASP truncates longer ones. Its
Psi4 reader, `interfaces/psi4/readfchk.py`, is Python 2: port its `print` statements to
Python 3 first.

Run:

    CAMCASP=/path/to/camcasp micromamba run -n cle-mbis python p2p_responses.py

Add molecule names (e.g. `formamide`) to run only those. Writes `p2p/<molecule>_000.p2p`
(with the CamCASP input and Psi4 log). About 30 min for all four on a laptop.

On macOS use the `CamCASP-6.1-macos` branch. Its binaries are x86_64 (they run under
Rosetta) and take their linear algebra from Apple's Accelerate, whose LU solver gives a
wrong density fit: formamide's response comes out 4% off. Supply the gcc runtime and
relink the binaries against OpenBLAS:

    cd $CAMCASP
    CONDA_SUBDIR=osx-64 micromamba create -y -p /tmp/x86 -c conda-forge libgfortran5 libopenblas=0.3.34
    mkdir runtime_libs
    cp -L /tmp/x86/lib/lib{gfortran.5,gomp.1,quadmath.0,gcc_s.1,gcc_s.1.1,omp,openblas.0}.dylib runtime_libs/
    for b in camcasp cluster process casimir pfit; do
      install_name_tool -change /System/Library/Frameworks/Accelerate.framework/Versions/A/Accelerate \
          /usr/local/opt/openblas/lib/libopenblas.0.dylib osx/gfortran/$b
      codesign -f -s - osx/gfortran/$b
      ln -s ../osx/gfortran/$b bin/$b
    done
    rm bin/no_psi4

# Polarizabilities

Drude polarizabilities fitted to the p2p responses (run `p2p_responses.py` first). Atoms of one
polarizable atom type share a polarizability; `SCHEME` at the top of the script picks the types:

- `"transferability"` (default): from the data. Every heavy atom starts as its own type, and two
  types of the same element are combined whenever that raises no molecule's fit error by more than
  `LAMBDA` (10%). The four monomers end up with five types for their 17 heavy atoms.
- `"transferability_conditioned"`: the same, but only atoms whose chemical environments match out
  to `RADIUS` bonds may combine.
- `"sapt2013"`: predefined types (one C and one N shared by imidazole and pyrazole, acetic acid's
  two O shared, every other heavy atom its own).

How the five types are chosen: two types are combined only if, after refitting, no molecule's
fit error (RMS over the grid) rises by more than 10%. Merging stops at five types because every
further merge would raise some molecule's error by at least 13.6%; any limit from 3.0% to 13.6%
gives the same five. Scoring all 5.5 million possible sharings confirmed that these five have the
lowest mean and worst-molecule error of any five-type model.

The fitting equations are in McDaniel's thesis (2014, Appendix A) and in McDaniel and Schmidt,
J. Phys. Chem. A 2013; the docstring of `polarizabilities.py` lists them.

Run:

    micromamba run -n cle-mbis python polarizabilities.py

Prints each type's polarizability and the atoms in it, and writes `polarizabilities.csv`.
About 25 min for `"transferability"`, 2 min for `"sapt2013"`.

# Force field XMLs

OpenMM Drude force fields from the MBIS charges and the polarizabilities (run
`partial_charges.py` and `polarizabilities.py` first). Run:

    micromamba run -n cle-mbis python force_field_xml.py

Writes `force_fields/<method>/<molecule>.xml` (one force field per charge method; MBIS is the
one used in this work) from the templates in `templates/`; only the charges and polarizabilities
change. Check them with:

    micromamba run -n cle-mbis python -m pytest test_force_field_xml.py

# Legacy Fortran fitter

The original Fortran fitter, kept for reference. Run:

    micromamba run -n cle-mbis bash legacy_fortran/run_fortran.sh formamide

Prints one polarizability per heavy atom. Needs gfortran.
