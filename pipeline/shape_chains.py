"""TIC 402710280.01 transit shape with the posterior chain kept, with limb darkening fixed (as in fit_new.py, which
also saves its own chain) or free (Kipping 2013 q1, q2, uniform priors): the sensitivity of Rp/R* and b to limb
darkening. Windows are centred exactly as in fit_new.py (audit_402710280.json, times_noS102_shape).

  python shape_chains.py fixed|free
Output: results/candidates/shape_chain_{fixed,free}.npy and shape_{fixed,free}.json
"""
import json
import sys

import env  # noqa: F401
import numpy as np

import transit_fit as tf
from audit_402710280 import beta_factor
from env import RESULTS

C = RESULTS / "candidates"


def main(mode):
    from fit_new import lcs_for
    au = json.load(open(RESULTS / "candidates" / "audit_402710280.json"))
    tcs = {t["sector"]: t["tc"] for t in au["times_noS102_shape"]}
    lcs = lcs_for(402710280)
    data = []
    for w in tf.windows_from_lcs(lcs, [tcs[s] for s in (39, 93, 100, 102)], 0.8, 5.3):
        oot = np.abs(w[0] - w[4]) > 0.75 * 5.3 / 24
        data.append((w[0], w[1], w[2] * beta_factor(w[0][oot], w[1][oot] - 1, w[2]), w[3], w[4]))
    out, ch = tf.fit_stacked_windows(data, 36.3399, 5182.0, nstep=30000, free_ld=(mode == "free"))
    np.save(C / f"shape_chain_{mode}.npy", ch.astype(np.float32))
    json.dump(out, open(C / f"shape_{mode}.json", "w"), indent=1, default=float)
    print(mode, {k: v for k, v in out.items() if k not in ("chain_summary",)}, out["chain_summary"]["steps_over_tau"])


if __name__ == "__main__":
    main(sys.argv[1])
