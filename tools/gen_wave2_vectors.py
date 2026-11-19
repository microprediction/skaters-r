"""Generate inst/parity/wave2_vectors.json from the Python reference.

Covers the public API added in the second porting wave: homogenize (every
coupling mode, custom candidates, over a spliced base), residual_transform
(m0, m1, m2) and ResidualDiagnostics, dantzig and search, PushforwardDist
affine/shift/scale, and the spec JSON and names. Same probes, burn-in and
input series as vectors.json.

Run from the R repository root with the reference source on the path:
    PYTHONPATH=/path/to/skaters/src python3 tools/gen_wave2_vectors.py

SKATERS_REF, if set, labels the reference in the output (for a checkout with
local changes).
"""
import json
import math
import os
import subprocess

import skaters
from skaters import (
    ResidualDiagnostics, dantzig, homogenize, homogenize_candidates, residual_transform, search,
    build, to_json, spec_name, ema_spec, ensemble_spec, conjugate_spec, diff_spec, frac_spec, std_spec,
)
from skaters.conjugate import conjugate
from skaters.dist import Dist
from skaters.leaf import leaf
from skaters.spec import leaf_spec, ema_t_spec
from skaters.tails import gpdtails
from skaters.transform import (
    ema_transform, gaussianize, log_transform, seasonal_scale, standardize, yeo_johnson,
)

V = json.load(open(os.path.join("inst", "parity", "vectors.json")))
BURN, PROBE, Q_LO, Q_HI = V["burn"], V["probe"], V["q_lo"], V["q_hi"]
MAIN = V["series"]
SERIES = {"main": MAIN, "positive": [math.exp(0.3 * y) for y in MAIN]}


def clean(x):
    if x != x:
        return "nan"
    if x == math.inf:
        return "inf"
    if x == -math.inf:
        return "-inf"
    return x


def probe(d):
    return [clean(v) for v in (d.mean, d.std, d.logpdf(PROBE), d.cdf(PROBE),
                               d.quantile(Q_LO), d.quantile(Q_HI), d.crps(PROBE))]


def spliced(k):
    return gpdtails(conjugate(leaf(k=k), ema_transform(0.1), k=k),
                    k=k, level=0.9, nexc=50, warmup=100)


def ema_leaf(k):
    return conjugate(leaf(k=k), ema_transform(0.1), k=k)


def scenarios():
    s = []
    s.append(("homogenize_leaf", "main", 1, homogenize(leaf(k=1))))
    s.append(("homogenize_ema", "main", 1, homogenize(ema_leaf(1))))
    for couple in ("none", "weights", "state", "shared"):
        s.append((f"homogenize_ema_k3_{couple}", "main", 3,
                  homogenize(ema_leaf(3), k=3, couple=couple, lam=0.4)))
    s.append(("homogenize_custom", "main", 1, homogenize(
        ema_leaf(1), candidates=homogenize_candidates((0.0, 0.9), (0.05,), (0.0, 0.3, 0.9), delta_max=0.6))))
    s.append(("homogenize_spliced", "main", 1, homogenize(spliced(1))))
    for model in ("m0", "m1", "m2"):
        s.append((f"residual_{model}", "main", 1, residual_transform(ema_leaf(1), model=model, lr=0.05)))
    s.append(("residual_m2_spliced", "main", 1, residual_transform(spliced(1), model="m2")))
    s.append(("dantzig", "main", 1, dantzig(k=1, expand_interval=40)))
    s.append(("search", "main", 1, search(k=1, expand_interval=60, max_pool=10)))
    # pushforward predictives handed to affine inverses
    s.append(("standardize_over_yj_exact", "main", 1, conjugate(
        conjugate(leaf(k=1), yeo_johnson(0.0, exact=True)), standardize(0.05))))
    s.append(("seasonal_scale_over_log", "positive", 3, conjugate(
        conjugate(leaf(k=3), log_transform(lower=-1.0), k=3), seasonal_scale(7), k=3)))
    s.append(("gaussianize_over_yj_exact", "main", 1, conjugate(
        conjugate(leaf(k=1), yeo_johnson(0.5, exact=True)), gaussianize())))
    return s


def affine_cases():
    inner = Dist([(0.6, 0.2, 0.7), (0.4, -0.5, 1.3)])
    _, inv_k = yeo_johnson(0.0, exact=True)
    _, inv_log = log_transform()
    bases = {"table": inv_k([inner], {})[0], "analytic": inv_log([inner], {})[0]}
    out = []
    for kind, d in bases.items():
        for a, b in ((2.5, -1.0), (0.3, 4.0), (-1.7, 0.5)):
            out.append({"kind": kind, "op": "affine", "a": a, "b": b, "probe": probe(d.affine(a, b))})
        out.append({"kind": kind, "op": "shift", "a": 1.5, "b": 0.0, "probe": probe(d.shift(1.5))})
        out.append({"kind": kind, "op": "scale", "a": -2.0, "b": 0.0, "probe": probe(d.scale(-2.0))})
        # an affine of an affine (reflected twice)
        out.append({"kind": kind, "op": "affine2", "a": -0.5, "b": 2.0,
                    "probe": probe(d.affine(-0.5, 2.0).affine(-3.0, 1.0))})
    return out


def diagnostics():
    rd = ResidualDiagnostics(K=4, halflife=20.0)
    zs = [(y - 0.5) / 2.0 for y in MAIN]
    rows = []
    for z in zs:
        rd.update(z)
        sm = rd.summary()
        rows.append([sm["n"], sm["E[z]"], sm["E[z2]"], sm["E[z3]"], sm["E[z4]"]] + sm["L"] + sm["cross_sq"])
    return {"K": 4, "halflife": 20.0, "z": zs, "rows": rows}


def specs():
    cases = [
        leaf_spec(1),
        leaf_spec(3),
        ema_spec(0.05, 1),
        ema_spec(1e-5, 1),
        ema_spec(1.0, 1),
        ema_spec(1, 1),
        ema_spec(0.1 + 0.2, 2),
        ema_spec(1e16, 1),
        ema_spec(123456.789, 1),
        ensemble_spec(ema_spec(0.01, 1), ema_spec(0.1, 1), k=1),
        conjugate_spec(ensemble_spec(ema_spec(0.01, 1), ema_spec(0.1, 1), k=1), diff_spec()),
        conjugate_spec(leaf_spec(1), frac_spec(0.4, 30)),
        conjugate_spec(leaf_spec(1), frac_spec(0.4)),
        conjugate_spec(leaf_spec(1), std_spec(0.05)),
        conjugate_spec(conjugate_spec(leaf_spec(2), ema_t_spec(1e-4)), diff_spec()),
        {**ensemble_spec(ema_spec(2.5e-7, 1), leaf_spec(1), k=1), "floor": 1e-8},
    ]
    out = []
    for c in cases:
        try:
            build(c)
            builds = True
        except (AssertionError, ValueError):
            builds = False                  # e.g. alpha >= 1 is no ema
        out.append({"json": to_json(c), "name": spec_name(c), "builds": builds})
    return out


def main():
    out = {"burn": BURN, "probe": PROBE, "q_lo": Q_LO, "q_hi": Q_HI, "series": SERIES, "scenarios": {}}
    src = os.path.dirname(os.path.dirname(os.path.dirname(skaters.__file__)))
    try:
        out["source"] = os.environ.get("SKATERS_REF") or subprocess.check_output(
            ["git", "-C", src, "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        out["source"] = "unknown"
    for name, series, k, sk in scenarios():
        state, rows = None, []
        for t, y in enumerate(SERIES[series]):
            dists, state = sk(y, state)
            if t >= BURN:
                rows.append([probe(d) for d in dists])
        out["scenarios"][name] = {"series": series, "k": k, "out": rows}
        print(f"[gen] {name}")
    out["affine"] = affine_cases()
    out["diagnostics"] = diagnostics()
    out["specs"] = specs()
    path = os.path.join("inst", "parity", "wave2_vectors.json")
    with open(path, "w") as f:
        json.dump(out, f, allow_nan=False)
    print(f"[gen] wrote {path} ({len(out['scenarios'])} scenarios, reference {out['source']})")


if __name__ == "__main__":
    main()
