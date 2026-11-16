#!/bin/bash
# Fit Drude polarizabilities with the original Fortran fitter (one type per heavy atom), for
# comparison with polarizabilities.py. Needs gfortran and p2p/<molecule>_000.p2p.
# Run from the tutorial root: bash legacy_fortran/run_fortran.sh [molecule ...]
# The fitter's inputs (<molecule>.clt, .constraint) are kept in legacy_fortran/inputs/.
set -e
SRC=$(cd "$(dirname "$0")" && pwd)
INPUTS=$SRC/inputs
WORK=$(mktemp -d)
mkdir -p "$INPUTS"
trap 'rm -rf "$WORK"' EXIT
# same compile line as SAPT-FF-CLE ff_fitting/pol/install/install_phoenix.sh:51
(cd "$WORK" && gfortran "$SRC"/variable.f90 "$SRC"/routines.f90 "$SRC"/scf_drude.f90 \
    "$SRC"/amoeba.f90 "$SRC"/main_p2p.f90 -o main_p2p_gfortran)
for mol in ${@:-formamide acetic_acid imidazole pyrazole}; do
  python "$SRC/fortran_inputs.py" "$mol" "$INPUTS"
  cp "$INPUTS/$mol.clt" "$INPUTS/$mol.constraint" "p2p/${mol}_000.p2p" "$WORK/"
  (cd "$WORK" && ./main_p2p_gfortran "$mol" > "$mol.out")
  # converged Drude charges q -> alpha = q^2 / k (k = 0.1 au), in bohr^3
  awk -v m="$mol" '/rms, max_d,min_d/{exit} f && NF==2 {printf "%-12s %-4s %10.4f bohr^3\n", m, $1, $2*$2/0.1}
                   /converged Drude charges/{f=1}' "$WORK/$mol.out"
done
