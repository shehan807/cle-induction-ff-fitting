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
