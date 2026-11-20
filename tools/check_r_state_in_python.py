"""Resume R-written checkpoints in the Python reference.

run_state_parity(out_dir = DIR) in the R port restores each reference
checkpoint, steps 20 ticks and writes its own state as DIR/<scenario>_<t>.json.
This script restores each of those in Python with state_from_json, steps the
next 20 ticks and compares the predictives with the reference's own (the
second half of `out` in the state vectors) at 1e-6.

    PYTHONPATH=/path/to/skaters/src:/path/to/skaters/parity \\
        python3 tools/check_r_state_in_python.py DIR
"""
import gzip
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_state_vectors as G  # noqa: E402
from skaters import state_from_json  # noqa: E402
from skaters.tails import SplicedDist  # noqa: E402

ATOL = RTOL = 1e-6


def num(e):
    return {"nan": math.nan, "inf": math.inf, "-inf": -math.inf}.get(e, e) if isinstance(e, str) else e


def close(got, want):
    if want != want:
        return got != got
    if math.isinf(want):
        return got == want
    return math.isfinite(got) and abs(got - want) <= ATOL + RTOL * abs(want)


def deep_tail(d, x):
    if not isinstance(d, SplicedDist):
        return False
    lower = sum(w * 0.5 * math.erfc(-(x - m) / (s * math.sqrt(2))) for w, m, s in d.body.components)
    upper = sum(w * 0.5 * math.erfc((x - m) / (s * math.sqrt(2))) for w, m, s in d.body.components)
    return min(lower, upper) < 1e-8


def main(out_dir):
    v = json.load(gzip.open(os.path.join("inst", "parity", "state_vectors.json.gz"), "rt"))
    cont = v["continue"]
    builders = {name: (k, sk, sname) for name, k, sk, sname in G.scenarios()}
    fails = checked = files = deep = 0
    for fn in sorted(os.listdir(out_dir)):
        if not fn.endswith(".json"):
            continue
        name, t = fn[:-5].rsplit("_", 1)
        t = int(t)
        k, sk, sname = builders[name]
        series = [num(e) for e in v["series"][sname]]
        cp = next(c for c in v["scenarios"][name]["checkpoints"] if c["t"] == t)
        state = state_from_json(json.load(open(os.path.join(out_dir, fn))))
        files += 1
        bad = 0
        for m in range(cont, 2 * cont):
            dists, state = sk(series[t + m], state)
            for h, d in enumerate(dists):
                got = G.probe(d)
                want = [num(e) for e in cp["out"][m][h]]
                dt = deep_tail(d, v["probe"])
                for j, (a, b) in enumerate(zip(got, want)):
                    checked += 1
                    a = num(a)
                    if dt and j in (2, 3) and not close(a, b) and math.isfinite(a) \
                            and abs(a - b) <= ATOL + 1e-4 * abs(b):
                        deep += 1
                        continue
                    if not close(a, b):
                        bad += 1
                        if bad <= 3:
                            print(f"FAIL {name} t={t} +{m + 1} h {h + 1} probe {j + 1}: got {a!r} want {b!r}")
        fails += bad
        print(("FAIL " if bad else "ok   ") + f"{name} t={t}")
    print(f"{files} R checkpoints resumed in Python, {checked} values checked, {deep} deep-tail probes at 1e-4")
    print("R->PYTHON RESUME FAILED" if fails or not files else "R->PYTHON RESUME OK")
    return 1 if fails or not files else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
