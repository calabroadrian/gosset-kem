"""
AFA-12 v12 — official lattice-estimator (malb/lattice-estimator, commit f31b86f) on key and ciphertext
LWE instances of ML-KEM and AFA-12 MaxSec. Mode 'rough' = core-SVP (0.292 beta); mode 'full' = default
cost model (MATZOV dual, dual-hybrid, primal hybrid, bdd), reported in log2 rop.
"""
import sys, math
sys.path.insert(0, 'lattice-estimator')
from estimator import LWE, ND
from estimator.lwe_parameters import LWEParameters

q = 3329
SETS = [("ML-KEM-512", 2, 3, 2, 10), ("AFA12-512-MaxSec", 2, 4, 3, 10),
        ("ML-KEM-768", 3, 2, 2, 10), ("AFA12-768-MaxSec", 3, 3, 2, 10),
        ("ML-KEM-1024", 4, 2, 2, 11), ("AFA12-1024-MaxSec", 4, 3, 2, 11)]

def instances(k, e1, e2, du):
    n = 256 * k
    sr = math.sqrt(e2 / 2 + (q / 2 ** du) ** 2 / 12)
    return {
        "key": LWEParameters(n=n, q=q, Xs=ND.CenteredBinomial(e1), Xe=ND.CenteredBinomial(e1), m=n),
        "ct": LWEParameters(n=n, q=q, Xs=ND.CenteredBinomial(e1), Xe=ND.CenteredBinomial(e2), m=n),
        "ct+round": LWEParameters(n=n, q=q, Xs=ND.CenteredBinomial(e1), Xe=ND.DiscreteGaussian(sr), m=n),
    }

mode = sys.argv[1] if len(sys.argv) > 1 else "rough"
for name, k, e1, e2, du in SETS:
    for inst, p in instances(k, e1, e2, du).items():
        if mode == "rough":
            r = LWE.estimate.rough(p, quiet=True)
        else:
            r = LWE.estimate(p, deny_list=("arora-gb", "bkw"), quiet=True)
        parts = {a: math.log2(v["rop"]) for a, v in r.items()}
        best = min(parts, key=parts.get)
        print(f"{name:18s} {inst:9s} min {parts[best]:6.1f} ({best})  " +
              "  ".join(f"{a}={v:.1f}" for a, v in sorted(parts.items())), flush=True)
