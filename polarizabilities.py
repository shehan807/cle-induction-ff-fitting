"""Fit Drude polarizabilities to the CamCASP p2p responses (run p2p_responses.py first).

Each heavy atom carries a Drude charge q on a spring k (alpha = q^2 / k). A probe charge at grid
point Q displaces the Drudes; the potential they create at grid point P is fitted to the CamCASP
response over all (P, Q) pairs. Atoms of one type share a polarizability. SCHEME picks the types:
  "transferability"              one type per heavy atom, then combine same-element types while
                                 no molecule's fit error rises by more than LAMBDA
  "transferability_conditioned"  the same, but only atoms whose chemical environments match out
                                 to RADIUS bonds may combine
  "sapt2013"                     one type per heavy atom, except acetic acid's two O and the
                                 imidazole + pyrazole C and N, which share

Equations (McDaniel PhD thesis 2014, Appendix A; molecule m, atoms a, grid points P, Q):
  alpha_a = q_a^2 / k                                                                   (A.14)
  model response  R~_PQ = sum_a q_a [1/|r_P - r_a - d_a(Q)| - 1/|r_P - r_a|]            (A.11)
  d_a(Q): Drude displacements where the force vanishes, k d_a = q_a E(r_a + d_a)   (min of A.12)
  Thole screening of Drude-Drude terms: 1 - (1 + s r/2) exp(-s r), s = p / (alpha_a alpha_b)^(1/6)
                                                                                        (A.13)
  fit error  RMS_m^2 = mean over P >= Q of (R_PQ - R~_PQ)^2                              (A.3)
  objective  sum over molecules of ln RMS_m
with k = 0.1 au and p = 2.0 (thesis p. 281; McDaniel and Schmidt, J. Phys. Chem. A 2013, 117, 2053).
The probe charge is -1, as in CamCASP (thesis p. 277).

How "transferability" picks the types: a merge is kept only if, after refitting,
  max over m of [ln RMS_m(after) - ln RMS_m(before)] < ln(1 + LAMBDA),
i.e. no molecule's error rises more than LAMBDA (10%) relative to before that merge. Merging stops
when no pair passes.

The Fortran fitter this replaces is in legacy_fortran/ ("Ref:" lines point to it).
Run: python polarizabilities.py
"""
import csv
import itertools
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.spatial.distance import cdist

from drude_fit.inputs import ANG_TO_BOHR, read_geometry

HERE = Path(__file__).resolve().parent
P2P = HERE / "p2p"
MOLECULES = ["formamide", "acetic_acid", "imidazole", "pyrazole"]   # same order as p2p_responses.py
SCHEME = "transferability"   # or "transferability_conditioned", "sapt2013"
LAMBDA = 0.10                # combine two types only if no molecule's fit error rises more than 10%
RADIUS = 1                   # transferability_conditioned: bonds out to which environments must match
COVALENT_RADII = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66}   # angstrom; Cordero et al. 2008
K = 0.1                      # Drude spring constant (au). Ref: legacy_fortran/variable.f90:14
THOLE = 2.0                  # Thole damping parameter. Ref: legacy_fortran/variable.f90:15
BOHR3_TO_A3 = ANG_TO_BOHR ** -3
START = {"C": 12.0, "N": 7.0, "O": 10.0}   # starting guesses (bohr^3), as in the production .constraint files


def read_p2p(path):
    """Grid (n, 3) in bohr and the response (n, n), lower triangle filled, from a .p2p file.

    Fortran unformatted records, each wrapped in 4-byte length markers (m0, m1): n, the n grid
    points, then the lower triangle row by row. Ref: legacy_fortran/routines.f90:62-96 (get_p2p)"""
    with open(path, "rb") as f:
        n = int(np.frombuffer(f.read(12), dtype="<i4")[1])
        grid = np.frombuffer(f.read(n * 32), dtype=[("m0", "<i4"), ("xyz", "<f8", 3), ("m1", "<i4")])
        tri = np.frombuffer(f.read(), dtype=[("m0", "<i4"), ("v", "<f8"), ("m1", "<i4")])
    response = np.zeros((n, n))
    response[np.tril_indices(n)] = tri["v"]
    return grid["xyz"].copy(), response


def thole(d, alpha_i, alpha_j):
    """Thole screening f(d) of the i-j pair (0 at contact, 1 far away) and df/dd.
    Thesis eq A.13. Ref: legacy_fortran/scf_drude.f90:215-268 (screen, dscreen)"""
    s = THOLE / (alpha_i * alpha_j) ** (1 / 6)
    e = np.exp(-s * d)
    return 1 - (1 + s * d / 2) * e, s * (1 + s * d / 2) * e - s / 2 * e


def screened_field(v, alpha_i, alpha_j):
    """Thole-screened field per unit source charge at separation v (Q, 3).
    Ref: legacy_fortran/scf_drude.f90:196-197, 203-204"""
    d = np.linalg.norm(v, axis=-1, keepdims=True)
    f, df = thole(d, alpha_i, alpha_j)
    return f * v / d ** 3 - df * v / d ** 2


def field_on_drudes(drudes, probes, cores, q):
    """Field at each Drude for every probe position, shape (Q, A, 3): from the probe (charge -1,
    since CamCASP's response is minus the potential a +1 charge induces) and from the other
    Drudes (q_j) and their cores (-q_j). Ref: legacy_fortran/scf_drude.f90:167-211 (field)"""
    v = drudes - probes[:, None, :]
    field = -v / np.linalg.norm(v, axis=-1, keepdims=True) ** 3
    alpha = q ** 2 / K
    for i in range(len(q)):
        for j in range(len(q)):
            if i != j:
                field[:, i] += q[j] * screened_field(drudes[:, i] - drudes[:, j], alpha[i], alpha[j])
                field[:, i] -= q[j] * screened_field(drudes[:, i] - cores[j], alpha[i], alpha[j])
    return field


def relax_drudes(probes, cores, q, tol=1e-8, max_steps=2000):
    """Drude positions where the field force q*E balances the spring k*(drude - core). Unsettled
    after max_steps (only at unphysically large trial alpha), the fit sees a large error there.
    Zero force of thesis eq A.12. Ref: legacy_fortran/scf_drude.f90:117-164 (findpositionshell,
    which mixes in 0.1 of each step)"""
    drudes = np.broadcast_to(cores, (len(probes), *cores.shape)).copy()
    for _ in range(max_steps):
        force = q[None, :, None] * field_on_drudes(drudes, probes, cores, q)
        if np.max(np.abs(force - K * (drudes - cores))) < tol:
            break
        drudes = cores + force / K
    return drudes


def rms_error(alpha, xyz, grid, response):
    """RMS difference between the Drude and CamCASP responses over the (P, Q <= P) pairs.
    Thesis eqs A.3, A.11. Ref: legacy_fortran/scf_drude.f90:5-39 (rms_p2p), 88-114
    (electrostatic_potential)"""
    on = alpha > 0                                   # hydrogens carry no Drude
    q = -np.sqrt(alpha[on] * K)                      # eq A.14. Ref: legacy_fortran/scf_drude.f90:52
    cores = xyz[on]
    drudes = relax_drudes(grid, cores, q)            # (Q, A, 3): one relaxed set per probe Q
    to_drude = cdist(grid, drudes.reshape(-1, 3)).reshape(len(grid), *drudes.shape[:2])  # (P, Q, A)
    to_core = cdist(grid, cores)                                                       # (P, A)
    drude_response = (q / to_drude).sum(-1) - (q / to_core).sum(-1)[:, None]
    rows, cols = np.tril_indices(len(grid))
    return np.sqrt(np.mean((drude_response[rows, cols] - response[rows, cols]) ** 2))


def load():
    """Labels, elements, coordinates (bohr), grid and response of each molecule."""
    return {name: (*read_geometry(HERE / "geometries" / f"{name}.pdb"),
                   *read_p2p(P2P / f"{name}_000.p2p")) for name in MOLECULES}


def ln_rms(name, alpha_of, data):
    """ln RMS error of one molecule (alpha_of: atom label -> alpha). Fits add these logs, so each
    molecule's error counts relative to its own size."""
    labels, elements, xyz, grid, response = data[name]
    alpha = np.array([alpha_of.get(lab, 0.0) for lab in labels])
    return np.log(rms_error(alpha, xyz, grid, response))


def fit(types, names, data, start, h=1e-4):
    """Fit the molecules `names` together: minimize their summed ln RMS over ln alpha (keeps alpha
    positive) with L-BFGS-B and finite-difference gradients. `types`: (molecule, label) -> type.
    Returns ({molecule: ln RMS}, {type: alpha}). Ref: legacy_fortran/amoeba.f90 (Nelder-Mead)"""
    names = sorted(names)                            # fixed order, so sums repeat exactly
    keys = sorted({t for (m, _), t in types.items() if m in names})
    type_of = {m: {lab: t for (mm, lab), t in types.items() if mm == m} for m in names}

    def loss(m, x):
        alpha = dict(zip(keys, np.exp(x)))
        return ln_rms(m, {lab: alpha[t] for lab, t in type_of[m].items()}, data)

    def objective(x):
        base = {m: loss(m, x) for m in names}
        grad = np.zeros(len(x))
        for i, t in enumerate(keys):
            step = x + h * (np.arange(len(x)) == i)
            grad[i] = sum(loss(m, step) - base[m] for m in names if t in type_of[m].values()) / h
        return sum(base.values()), grad

    res = minimize(objective, np.log([start[t] for t in keys]), jac=True, method="L-BFGS-B")
    # finite-difference noise can stop the line search at the minimum without "success"
    assert res.success or np.abs(res.jac).max() < 1e-3, f"{names}: {res.message}"
    return {m: loss(m, res.x) for m in names}, dict(zip(keys, np.exp(res.x)))


def linked(types, names):
    """`names` plus every molecule tied to them through a shared atom type."""
    names = set(names)
    while True:
        shared = {t for (m, _), t in types.items() if m in names}
        more = {m for (m, _), t in types.items() if t in shared} - names
        if not more:
            return sorted(names)
        names |= more


def fit_types(types, data):
    """Fit fixed atom types; molecules that share a type are fitted together.
    Returns ({type: alpha}, {molecule: ln RMS})."""
    alpha, now, left = {}, {}, set(data)
    while left:
        names = linked(types, {left.pop()})
        left -= set(names)
        loss, a = fit(types, names, data, {t: START[t[0]] for t in types.values()})
        alpha.update(a)
        now.update(loss)
    return alpha, now


def sapt2013_types(data):
    """The predefined types (see the module docstring)."""
    types = {}
    for m, (labels, elements, *_) in data.items():
        for lab, el in zip(labels, elements):
            if el == "H":
                continue
            if m in ("imidazole", "pyrazole"):
                types[(m, lab)] = f"{el}: imidazole + pyrazole"
            elif m == "acetic_acid" and el == "O":
                types[(m, lab)] = "O: acetic_acid"
            else:
                types[(m, lab)] = f"{el}: {m} {lab}"
    return types


def environment(data, radius):
    """Morgan-style label per heavy atom: element and number of bonded atoms (radius 0), then the
    sorted labels of its bonded atoms, once per extra bond of radius. Bonded = closer than 1.2x the
    summed covalent radii."""
    env_of = {}
    for m, (labels, elements, xyz, *_) in data.items():
        r = np.array([COVALENT_RADII[e] for e in elements]) * ANG_TO_BOHR
        bonded = (cdist(xyz, xyz) < 1.2 * (r[:, None] + r)) & ~np.eye(len(elements), dtype=bool)
        env = [f"{e}{n}" for e, n in zip(elements, bonded.sum(1))]
        for _ in range(radius):
            env = [env[i] + "(" + ",".join(sorted(env[j] for j in np.flatnonzero(bonded[i]))) + ")"
                   for i in range(len(env))]
        env_of.update({(m, lab): e for lab, e, el in zip(labels, env, elements) if el != "H"})
    return env_of


def curvature(name, alpha_of, data, h=0.05):
    """Quadratic model of one molecule's ln RMS around its own best fit, in ln alpha:
    (labels, ln alpha at the minimum, ln RMS there, Hessian by central differences)."""
    labels = list(alpha_of)
    x0 = np.log([alpha_of[lab] for lab in labels])

    def f(x):
        return ln_rms(name, dict(zip(labels, np.exp(x))), data)

    n, e, f0 = len(x0), h * np.eye(len(x0)), f(x0)
    hess = np.zeros((n, n))
    for i in range(n):
        hess[i, i] = (f(x0 + e[i]) - 2 * f0 + f(x0 - e[i])) / h ** 2
        for j in range(i):
            hess[i, j] = hess[j, i] = (f(x0 + e[i] + e[j]) - f(x0 + e[i] - e[j])
                                       - f(x0 - e[i] + e[j]) + f(x0 - e[i] - e[j])) / (4 * h ** 2)
    return labels, x0, f0, hess


def predict(types, names, model):
    """Summed ln RMS of `names` under `types` from the quadratic models, and the best alpha per
    type. The constrained minimum of a sum of quadratics is one linear solve."""
    names = sorted(names)
    keys = sorted({t for (m, _), t in types.items() if m in names})
    a, b, maps = np.zeros((len(keys), len(keys))), np.zeros(len(keys)), {}
    for m in names:
        labels, x0, _, hess = model[m]
        maps[m] = np.array([[types[(m, lab)] == t for t in keys] for lab in labels], dtype=float)
        a += maps[m].T @ hess @ maps[m]
        b += maps[m].T @ hess @ x0
    y = np.linalg.solve(a, b)
    total = 0.0
    for m in names:
        labels, x0, f0, hess = model[m]
        d = maps[m] @ y - x0
        total += f0 + d @ hess @ d / 2
    return total, dict(zip(keys, np.exp(y)))


def combine(types, t1, t2):
    """`types` with t1 and t2 made one type, and its name (element and the atoms in it)."""
    members = sorted(atom for atom, t in types.items() if t in (t1, t2))
    new = t1.split(":")[0] + ": " + ", ".join(f"{m} {lab}" for m, lab in members)
    return {atom: new if t in (t1, t2) else t for atom, t in types.items()}, new


def transferability_types(data, radius=None):
    """Start from one type per heavy atom. Each round, rank the same-element pairs of types by
    their predicted error rise (curvature model), refit them exactly in that order and keep the
    first that raises no molecule's error by more than LAMBDA. Stop when none passes. With `radius`,
    only types whose environments match (see environment) may combine.
    Returns (types, {type: alpha}, {molecule: ln RMS})."""
    types = {(m, lab): f"{el}: {m} {lab}" for m, (labels, elements, *_) in data.items()
             for lab, el in zip(labels, elements) if el != "H"}
    alpha, now = fit_types(types, data)
    model = {m: curvature(m, {lab: alpha[t] for (mm, lab), t in types.items() if mm == m}, data)
             for m in data}
    limit = np.log1p(LAMBDA)                         # LAMBDA as a rise in ln RMS
    env_of = environment(data, radius) if radius is not None else None
    tried = set()
    while True:
        env = {t: env_of[a] for a, t in types.items()} if env_of else {}   # one label per type
        options = []
        for t1, t2 in itertools.combinations(sorted(set(types.values())), 2):
            if (t1.split(":")[0] != t2.split(":")[0] or env.get(t1) != env.get(t2)
                    or (t1, t2) in tried):
                continue
            trial, new = combine(types, t1, t2)
            gaining = {m for (m, _), t in types.items() if t in (t1, t2)}
            names = linked(trial, gaining)
            after, start = predict(trial, names, model)
            before, _ = predict(types, names, model)
            options.append(((after - before) / len(gaining), t1, t2, new, trial, names, start))

        combined = False
        for predicted, t1, t2, new, trial, names, start in sorted(options, key=lambda o: o[0]):
            if predicted > 1.5 * limit:              # predictions are good to about limit / 2
                break
            tried.add((t1, t2))
            loss, a = fit(trial, names, data, start)
            rise = max(loss[m] - now[m] for m in names)
            print(f"{'combined  ' if rise < limit else 'kept apart'}  error {100 * np.expm1(rise):+6.2f}%  {new}",
                  flush=True)
            if rise < limit:
                types = trial
                now.update(loss)
                alpha = {t: v for t, v in alpha.items() if t in types.values()}
                alpha.update(a)
                combined = True
                break
        if not combined:
            return types, alpha, now


def write2csv(types, alpha, path="polarizabilities.csv"):
    """Write (molecule, label, atom type, alpha in bohr^3 and A^3) rows to `path`."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(("molecule", "label", "atom_type", "alpha_bohr3", "alpha_A3"))
        w.writerows((m, lab, t, alpha[t], alpha[t] * BOHR3_TO_A3) for (m, lab), t in types.items())


def main():
    data = load()
    if SCHEME == "sapt2013":
        types = sapt2013_types(data)
        alpha, now = fit_types(types, data)
    else:
        radius = RADIUS if SCHEME == "transferability_conditioned" else None
        types, alpha, now = transferability_types(data, radius)
    print(f"\n{len(alpha)} polarizable atom types ({SCHEME})")
    for t, a in sorted(alpha.items()):
        print(f"  {a:8.4f} bohr^3  {a * BOHR3_TO_A3:.4f} A^3   {t}")
    for m in MOLECULES:
        print(f"  {m:12s} RMS {np.exp(now[m]):.4e}")
    write2csv(types, alpha)


if __name__ == "__main__":
    main()
