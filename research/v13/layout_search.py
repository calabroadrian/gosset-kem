"""
AFA-12 v13 — search a public coefficient layout for the 32 E8 blocks such that, inside every block, the 28
pairwise lags are all distinct (lag class = min(d, 256-d); negacyclic: lag 256-d is -A(d)).
Goal: an attacker who grinds ciphertexts for a large autocorrelation at one lag can only raise ONE pair per
block (weight 1/4 in the type-B directions), instead of up to 7 pairs with the consecutive layout.
Secondary goal: spread each lag class over as few blocks as possible (balance).
"""
import numpy as np, itertools, json, sys

N = 256
def lagc(a, b):
    d = (b - a) % N
    return min(d, N - d)

def search(seed):
    rng = np.random.default_rng(seed)
    free = list(range(N)); rng.shuffle(free); free = set(free)
    use = np.zeros(129, int)                 # how many blocks use each lag class
    blocks = []
    for b in range(32):
        best = None
        for attempt in range(4000):
            pool = sorted(free); start = pool[0] if attempt == 0 else rng.choice(pool)
            blk = [start]; lags = set()
            cand = [x for x in pool if x != start]; rng.shuffle(cand)
            # greedy: add positions keeping lags distinct, prefer rarely-used lag classes
            cand.sort(key=lambda x: use[lagc(start, x)] + rng.random())
            for x in cand:
                nl = [lagc(y, x) for y in blk]
                if len(set(nl)) < len(nl) or any(l in lags for l in nl): continue
                blk.append(x); lags.update(nl)
                if len(blk) == 8: break
            if len(blk) == 8:
                cost = sum(use[l] for l in lags)
                if best is None or cost < best[0]: best = (cost, blk, lags)
                if attempt > 200: break
        if best is None: return None
        _, blk, lags = best
        for x in blk: free.discard(x)
        for l in lags: use[l] += 1
        blocks.append(sorted(int(x) for x in blk))
    return blocks, use

if __name__ == "__main__":
    for seed in range(200):
        r = search(seed)
        if r is None: continue
        blocks, use = r
        print(f"seed {seed}: found. lag-class usage max {use[1:].max()}, mean {use[1:].mean():.2f}, lag 1 used by {use[1]} blocks")
        json.dump({"blocks": blocks}, open("layout.json", "w"))
        break
    else:
        print("no layout found"); sys.exit(1)
    # verify
    allpos = sorted(x for b in blocks for x in b); assert allpos == list(range(N))
    for b in blocks:
        ls = [lagc(x, y) for x, y in itertools.combinations(b, 2)]; assert len(set(ls)) == 28
    print("verified: partition of 0..255, 28 distinct lag classes in every block")
    print("first blocks:", blocks[:3])
