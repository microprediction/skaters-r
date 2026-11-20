"""Generate inst/parity/transform_vectors.json from the Python reference.

Covers the transforms and pushforward paths that parity/vectors.json in the
reference does not: seasonal_scale, gaussianize (all guards, smooth),
log_transform, flesaker, the exact flag on power_transform and yeo_johnson,
a pushforward fed through another component-reading transform, and every
component-reading transform over a tail-spliced predictive. Same probes,
burn-in and input series as vectors.json.

Run from the R repository root with the reference source on the path:
    PYTHONPATH=/path/to/skaters/src python3 tools/gen_transform_vectors.py
"""
import json
import math
import os
import subprocess

import skaters
from skaters.api import flesaker
from skaters.conjugate import conjugate
from skaters.leaf import leaf
from skaters.tails import gpdtails
from skaters.transform import (
    ema_transform, gaussianize, log_transform, power_transform, seasonal_anchor,
    seasonal_difference, seasonal_scale, yeo_johnson,
)

V = json.load(open(os.path.join("inst", "parity", "vectors.json")))
BURN, PROBE, Q_LO, Q_HI = V["burn"], V["probe"], V["q_lo"], V["q_hi"]
MAIN = V["series"]
SERIES = {
    "main": MAIN,
    "positive": [math.exp(0.3 * y) for y in MAIN],
}


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
    # GPD splice live from step 100 of the 150-point series.
    return gpdtails(conjugate(leaf(k=k), ema_transform(0.1), k=k),
                    k=k, level=0.9, nexc=50, warmup=100)


# Scenarios whose reference values are float-limited, with the relative
# tolerance they are checked at. A spliced predictive's quantile past the
# upper splice is the body's quantile at 1 - 1e-12, where the body CDF is
# flat to within its last bit, so the solved z is determined only to about
# 1e-6 absolute. Exponential amplification (yeo_johnson at lambda 0) or a
# table window placed by that quantile (gaussianize) carries it to the probes.
RTOL = {"spliced_yj_log": 1e-5, "spliced_yj_log_k3": 1e-5,
        "spliced_gaussianize": 1e-5, "spliced_gaussianize_k3": 1e-5}


def scenarios():
    s = []
    for k in (1, 3):
        suf = "" if k == 1 else f"_k{k}"

        def add(name, series, sk):
            s.append((name + suf, series, k, sk))

        add("seasonal_scale", "main", conjugate(leaf(k=k), seasonal_scale(7), k=k))
        add("seasonal_scale_center", "main", conjugate(leaf(k=k), seasonal_scale(7, center=True), k=k))
        add("seasonal_scale_sdiff", "main", conjugate(
            conjugate(leaf(k=k), seasonal_difference(7), k=k), seasonal_scale(7), k=k))
        add("gaussianize", "main", conjugate(leaf(k=k), gaussianize(), k=k))
        add("gaussianize_smooth", "main", conjugate(leaf(k=k), gaussianize(smooth=True), k=k))
        add("gaussianize_delete", "main", conjugate(leaf(k=k), gaussianize(guard="delete", slope_cap=2.0), k=k))
        add("gaussianize_emin", "main", conjugate(leaf(k=k), gaussianize(guard="emin", slope_cap=2.0), k=k))
        add("gaussianize_move_tight", "main", conjugate(leaf(k=k), gaussianize(slope_cap=2.0, max_knots=15), k=k))
        add("power_exact", "main", conjugate(leaf(k=k), power_transform(0.5, exact=True), k=k))
        add("yeojohnson_log_exact", "main", conjugate(leaf(k=k), yeo_johnson(0.0, exact=True), k=k))
        add("yeojohnson_half_exact", "main", conjugate(leaf(k=k), yeo_johnson(0.5, exact=True), k=k))
        # a pushforward read through the components view by the next transform
        add("power_over_yj_exact", "main", conjugate(
            conjugate(leaf(k=k), yeo_johnson(0.5, exact=True), k=k), power_transform(0.5), k=k))
        add("yj_over_power_exact", "main", conjugate(
            conjugate(leaf(k=k), power_transform(0.5, exact=True), k=k), yeo_johnson(0.0), k=k))
        add("log_lower", "positive", conjugate(leaf(k=k), log_transform(), k=k))
        add("log_lower_shifted", "main", conjugate(leaf(k=k), log_transform(lower=-6.0), k=k))
        add("log_upper", "main", conjugate(leaf(k=k), log_transform(upper=10.0), k=k))
        # tail-spliced predictives through every component-reading transform
        add("spliced_power", "main", conjugate(spliced(k), power_transform(0.5), k=k))
        add("spliced_yj_log", "main", conjugate(spliced(k), yeo_johnson(0.0), k=k))
        add("spliced_yj_half", "main", conjugate(spliced(k), yeo_johnson(0.5), k=k))
        add("spliced_gaussianize", "main", conjugate(spliced(k), gaussianize(), k=k))
        add("spliced_ema", "main", conjugate(spliced(k), ema_transform(0.05), k=k))
        add("spliced_seasonal_difference", "main", conjugate(spliced(k), seasonal_difference(2), k=k))
        add("spliced_seasonal_anchor", "main", conjugate(spliced(k), seasonal_anchor(2), k=k))
        add("spliced_seasonal_scale", "main", conjugate(spliced(k), seasonal_scale(7), k=k))
        add("spliced_log", "positive", conjugate(spliced(k), log_transform(), k=k))
        add("spliced_power_over_yj", "main", conjugate(
            conjugate(spliced(k), yeo_johnson(0.5), k=k), power_transform(0.5), k=k))
    s.append(("flesaker", "positive", 1, flesaker(k=1)))
    return s


def main():
    out = {"burn": BURN, "probe": PROBE, "q_lo": Q_LO, "q_hi": Q_HI,
           "series": SERIES, "scenarios": {}}
    src = os.path.dirname(os.path.dirname(os.path.dirname(skaters.__file__)))
    try:
        out["source"] = subprocess.check_output(
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
        if name in RTOL:
            out["scenarios"][name]["rtol"] = RTOL[name]
        print(f"[gen] {name}")
    path = os.path.join("inst", "parity", "transform_vectors.json")
    with open(path, "w") as f:
        json.dump(out, f, allow_nan=False)
    print(f"[gen] wrote {path} ({len(out['scenarios'])} scenarios, reference {out['source']})")


if __name__ == "__main__":
    main()
