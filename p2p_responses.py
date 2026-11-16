"""Compute the CamCASP point-to-point (p2p) response for four monomers:
acetic acid, formamide, imidazole, pyrazole. Geometries in geometries/, GRAC shifts in
grac_shifts.csv (run grac_shifts.py first).

The response is the change in electrostatic potential at each of 2000 grid points around
the molecule when a charge sits at another grid point (PBE0/d-aug-cc-pVTZ, ALDA+CHF
kernel, GRAC asymptotic correction). polarizabilities.py fits Drude oscillators to it.

Needs CamCASP 6.1 with $CAMCASP set to its install (on macOS relinked to OpenBLAS; see README).
Env: environment.yml.  Run: python p2p_responses.py [molecule ...]   (~30 min for all four)
"""
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

from drude_fit.inputs import read_geometry

HERE = Path(__file__).resolve().parent
CAMCASP = Path(os.environ["CAMCASP"])
MOLECULES = ["formamide", "acetic_acid", "imidazole", "pyrazole"]   # cheapest first
Z = {"H": 1, "C": 6, "N": 7, "O": 8}
MEMORY_GB = 24

# CamCASP input (.clt); davtz = d-aug-cc-pVTZ
CLT = """Title {name} : davtz properties
!
Global
  Units Bohr
  Overwrite Yes
End

Molecule {name}
I.P. {ip:.8f} au
HOMO {homo:.8f} au
Charge 0
{atoms}
End

Run
  Properties
  Molecule {name}
  Basis davtz
  File-prefix {name}
  Functional PBE0
  Kernel ALDA+CHF
  Orient files
  Process file
  Interface files
  Sites file
  Memory {memory} GB
End

Finish
"""


def write_clt(name, run_dir):
    """CamCASP input: geometry (bohr), the GRAC IP and HOMO, and the level of theory."""
    labels, elements, xyz = read_geometry(HERE / "geometries" / f"{name}.pdb")
    with open(HERE / "grac_shifts.csv") as f:
        grac = next(row for row in csv.DictReader(f) if row["molecule"] == name)
    atoms = "\n".join(f"{lab:<4} {Z[el]:>6.1f} {x:>16.8f} {y:>16.8f} {z:>16.8f}   Type {el}"
                      for lab, el, (x, y, z) in zip(labels, elements, xyz))
    clt = run_dir / f"{name}.clt"
    clt.write_text(CLT.format(name=name, ip=float(grac["ip_au"]), homo=float(grac["homo_au"]),
                              atoms=atoms, memory=MEMORY_GB))
    return clt


def camcasp_env(scratch):
    """Environment CamCASP needs: its own paths, the Psi4 in this env, and scratch space."""
    env = dict(os.environ)
    env.update(
        PATH=f"{CAMCASP / 'bin'}:{Path(sys.executable).parent}:{env['PATH']}",
        ARCH="osx" if sys.platform == "darwin" else "x86-64",
        PSI4_HOME=sys.prefix + "/",
        PSIPATH=f"{CAMCASP / 'basis/psi4'}:{CAMCASP / 'basis/psi4/for-psi4-lib'}",
        SCRATCH=str(scratch), PSI_SCRATCH=str(scratch),
        DYLD_FALLBACK_LIBRARY_PATH=str(CAMCASP / "runtime_libs"),   # macOS: gfortran + OpenBLAS
    )
    return env


def run_camcasp(name, run_dir, env):
    """Run CamCASP's `cluster` setup, the Psi4 SCF, then `camcasp`; return the .p2p file
    holding the static and 10 imaginary-frequency responses."""
    clt = write_clt(name, run_dir)
    subprocess.run([sys.executable, CAMCASP / "bin/runcamcasp.py", name, "--clt", clt,
                    "-d", run_dir / name, "--work", name, "--scfcode", "psi4",
                    "-M", str(MEMORY_GB), "--ifexists", "delete"], cwd=run_dir, env=env)
    found = list((run_dir / name / "OUT").glob("*_lim*_f*.p2p"))
    if len(found) != 1:
        raise RuntimeError(f"{name}: expected one .p2p from CamCASP, found {len(found)}; "
                           f"see the logs in {run_dir / name / 'OUT'}")
    return found[0]


def split_static(name, combined, env):
    """Split the combined file into one file per frequency; return the static one (_000)."""
    n_imaginary = int(combined.stem.rsplit("_f", 1)[1]) - 1   # name ends _f<n>: static + imaginary
    process_input = (f"read P2P pols for {name}\n    Frequencies Static + {n_imaginary}\n"
                     f"    P2P-Pols {combined.name} SPLIT\n  End\n")
    with open(combined.parent / "split_p2p.out", "w") as log:
        subprocess.run(["process"], input=process_input, text=True, stdout=log, cwd=combined.parent,
                       env=env, check=True)
    return combined.parent / f"{name}_000.p2p"


def main():
    out = HERE / "p2p"
    out.mkdir(exist_ok=True)
    for name in sys.argv[1:] or MOLECULES:
        run_dir = HERE / "camcasp_work" / name
        shutil.rmtree(run_dir, ignore_errors=True)
        (run_dir / "scratch").mkdir(parents=True)
        env = camcasp_env(run_dir / "scratch")
        static = split_static(name, run_camcasp(name, run_dir, env), env)
        shutil.copy(static, out)
        shutil.copy(run_dir / f"{name}.clt", out)                   # keep the input and the SCF log
        shutil.copy(run_dir / name / "OUT" / f"{name}_A.out", out)
        shutil.rmtree(run_dir)                                      # ~10 GB of intermediates
        print(f"{name:12s} -> p2p/{static.name}")


if __name__ == "__main__":
    main()
