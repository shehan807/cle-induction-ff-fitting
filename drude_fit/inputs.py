"""Turn a molecule's PDB into the arrays the fit needs."""

from pathlib import Path

import numpy as np

ANG_TO_BOHR = 1.8897259886  # coordinates come in Angstrom; the physics works in Bohr


def read_geometry(pdb_path):
    """Read a PDB into (labels, elements, xyz_bohr).

    labels and elements are one per atom; xyz is (n_atoms, 3) in Bohr, rounded to 8
    decimals so the atoms sit exactly where the quantum response was computed.
    """
    labels, elements, xyz = [], [], []
    for line in Path(pdb_path).read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        col = line.split()
        labels.append(col[2])
        elements.append(col[-1])
        xyz.append([round(float(col[c]) * ANG_TO_BOHR, 8) for c in (5, 6, 7)])
    return labels, elements, np.array(xyz)
