"""Generate inst/parity/state_vectors.json.gz: cross-language checkpoints.

For every exported skater the R port carries (the scenarios of the
reference's tests/test_state_json.py, less mahalanobis, zbank and gpdtail,
slated for extraction, and conforming_leaf and ensemble_leaf, which the
reference does not export), Python runs the series and writes, at
each checkpoint t (60, 250, 530 where the series allows):

    state        the canonical JSON state after t observations
    out          probes of the next 40 predictives (steps t+1 .. t+40)
    state_next   the canonical JSON state after t+20 observations

R restores `state`, steps 20 and must match out[0:20] to 1e-6; its own
state then is compared with `state_next`, written as JSON, and handed back to
Python (tools/check_r_state_in_python.py), which steps it 20 more and must
match out[20:40].

Run from the R repository root with the reference source on the path:
    PYTHONPATH=/path/to/skaters/src:/path/to/skaters/parity python3 tools/gen_state_vectors.py
SKATERS_REF, if set, labels the reference in the output; STATE_CHECKPOINTS
(comma-separated) and STATE_VECTORS_OUT override the checkpoints and path.
"""
import gzip
import json
import math
import os
import subprocess

import gen_vectors  # the reference's parity/gen_vectors.py

import skaters
from skaters import (
    laplace, granger, flesaker, dantzig, multiscale, parade, sticky, homogenize,
    residual_transform, conjugate, leaf, ema, precision_weighted_ensemble,
    bayesian_ensemble, terminal_leaf_ensemble, search, state_to_json,
)
from skaters import transform as T
from skaters.leaf import scale_mixture_leaf, crps_leaf, garch_leaf
from skaters.tails import gpdtails

N = 620
CHECKPOINTS = tuple(int(c) for c in os.environ.get("STATE_CHECKPOINTS", "60,250,530").split(","))
CONTINUE = 20
PROBE, Q_LO, Q_HI = 0.3, 0.1, 0.9


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


def _c(t, k=1):
    return conjugate(leaf(k=k), t, k=k)


def every_transform():
    return [
        ("difference", T.difference()),
        ("fractional_difference", T.fractional_difference(0.4, 30)),
        ("standardize", T.standardize()),
        ("ema_transform", T.ema_transform(0.1)),
        ("ou_transform", T.ou_transform(0.1)),
        ("theta", T.theta(0.1)),
        ("drift", T.drift(0.05, 0.01)),
        ("holt_linear", T.holt_linear(0.1, 0.05)),
        ("garch", T.garch()),
        ("seasonal_difference", T.seasonal_difference(7)),
        ("seasonal_anchor", T.seasonal_anchor(7)),
        ("seasonal_scale", T.seasonal_scale(7)),
        ("seasonal_scale_center", T.seasonal_scale(7, center=True)),
        ("power_transform", T.power_transform(0.5)),
        ("power_transform_exact", T.power_transform(0.5, exact=True)),
        ("yeo_johnson", T.yeo_johnson(0.5)),
        ("yeo_johnson_exact", T.yeo_johnson(0.5, exact=True)),
        ("ar", T.ar(2)),
        ("grouped_ar", T.grouped_ar(8)),
        ("gaussianize", T.gaussianize()),
        ("gaussianize_smooth", T.gaussianize(smooth=True)),
    ]


SERIES = {
    "main": gen_vectors.make_series(n=N, seed=12345),
    "repeat": gen_vectors.make_repeat_series(n=N, seed=99),
    "gaps": gen_vectors.make_gap_series(n=N, seed=4242),
}
SERIES["short"] = SERIES["main"][:300]
SERIES["positive"] = [abs(v) + 1.0 for v in SERIES["main"]]


def scenarios():
    s = []
    for name, t in every_transform():
        s.append(("conjugate_" + name, 1, _c(t), "positive" if name.startswith("power") else "main"))
    s.append(("conjugate_log_transform", 1, conjugate(leaf(k=1), T.log_transform(), k=1), "positive"))
    s.append(("laplace", 1, laplace(1), "main"))
    s.append(("laplace_k3", 3, laplace(3), "main"))
    s.append(("laplace_likelihood", 1, laplace(1, objective="likelihood"), "main"))
    s.append(("laplace_gaussian_tails", 1, laplace(1, tails="gaussian"), "short"))
    s.append(("laplace_repeat", 1, laplace(1), "repeat"))
    s.append(("laplace_gaps", 1, laplace(1), "gaps"))
    s.append(("laplace_k3_gaps", 3, laplace(3), "gaps"))
    s.append(("granger", 1, granger(1), "main"))
    s.append(("flesaker", 1, flesaker(1, lower=-100.0), "main"))
    s.append(("dantzig", 1, dantzig(1, expand_interval=50), "short"))
    s.append(("search_k2", 2, search(k=2, expand_interval=40, max_pool=12), "short"))
    s.append(("multiscale_k4", 4, multiscale(lambda kk: _c(T.ema_transform(0.1), kk), 4), "short"))
    s.append(("parade_spliced", 1, parade(gpdtails(_c(T.ema_transform(0.1)), k=1, nexc=50, warmup=100), 1), "main"))
    s.append(("parade_laplace_k3", 3, parade(laplace(3), 3), "main"))
    s.append(("sticky_ema_repeat", 1, sticky(_c(T.ema_transform(0.1)), k=1), "repeat"))
    s.append(("homogenize", 1, homogenize(_c(T.ema_transform(0.1)), k=1), "short"))
    s.append(("homogenize_k2_weights", 2, homogenize(_c(T.ema_transform(0.1), 2), k=2, couple="weights"), "short"))
    s.append(("homogenize_k2_shared", 2, homogenize(_c(T.ema_transform(0.1), 2), k=2, couple="shared"), "short"))
    for model in ("m0", "m1", "m2"):
        s.append((f"residual_transform_{model}", 1, residual_transform(_c(T.ema_transform(0.1)), model=model), "short"))
    s.append(("pw_ensemble_k2", 2, precision_weighted_ensemble([ema(0.05, k=2), ema(0.2, k=2)], k=2), "short"))
    s.append(("bayes_ensemble_k2", 2, bayesian_ensemble([ema(0.05, k=2), _c(T.difference(), 2)], k=2), "short"))
    s.append(("terminal_adaptive", 1, terminal_leaf_ensemble(
        [ema(0.05, k=1), _c(T.difference())], k=1, depths=[0, 1], adaptive_temperature=True), "short"))
    s.append(("terminal_k3", 3, terminal_leaf_ensemble(
        [ema(0.05, k=3), _c(T.difference(), 3)], k=3, depths=[0, 1], forget=0.99), "short"))
    for nm, lf in (("leaf", leaf), ("scale_mixture_leaf", scale_mixture_leaf), ("crps_leaf", crps_leaf),
                   ("garch_leaf", garch_leaf)):
        s.append((nm, 1, lf(k=1), "main"))
    s.append(("ema_k2", 2, ema(0.1, k=2), "short"))
    return s


def main():
    out = {"probe": PROBE, "q_lo": Q_LO, "q_hi": Q_HI, "continue": CONTINUE,
           "series": {k: [clean(v) for v in v_] for k, v_ in SERIES.items()}, "scenarios": {}}
    src = os.path.dirname(os.path.dirname(os.path.dirname(skaters.__file__)))
    try:
        out["source"] = os.environ.get("SKATERS_REF") or subprocess.check_output(
            ["git", "-C", src, "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        out["source"] = "unknown"
    for name, k, sk, sname in scenarios():
        series = SERIES[sname]
        n = len(series)
        cps = [c for c in CHECKPOINTS if c + 2 * CONTINUE < n]
        state, rows, cp_out = None, [], []
        outs = {}
        states = {}
        for t in range(n):
            if t in cps:
                states[t] = json.loads(json.dumps(state_to_json(state), allow_nan=False))
            if t - CONTINUE in cps:
                states[(t - CONTINUE, "next")] = json.loads(json.dumps(state_to_json(state), allow_nan=False))
            dists, state = sk(series[t], state)
            for c in cps:
                if c <= t < c + 2 * CONTINUE:
                    outs.setdefault(c, []).append([probe(d) for d in dists])
            if cps and t >= max(cps) + 2 * CONTINUE:
                break
        for c in cps:
            cp_out.append({"t": c, "state": states[c], "state_next": states[(c, "next")], "out": outs[c]})
        out["scenarios"][name] = {"series": sname, "k": k, "checkpoints": cp_out}
        print(f"[gen] {name}: checkpoints {cps}")
    path = os.environ.get("STATE_VECTORS_OUT") or os.path.join("inst", "parity", "state_vectors.json.gz")
    with gzip.open(path, "wt") as f:
        json.dump(out, f, allow_nan=False)
    print(f"[gen] wrote {path} ({len(out['scenarios'])} scenarios, reference {out['source']})")


if __name__ == "__main__":
    main()
