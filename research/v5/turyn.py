"""Find H' (permuted Hamming) so that Turyn's construction gives the Golay code."""
import itertools, numpy as np
G = np.array([[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[0,1,0,1,0,1,0,1]])
def words(G): return [int(''.join(map(str,(np.array(b)@G)%2)),2) for b in itertools.product([0,1],repeat=4)]
H = words(G)
def perm_word(w, p): bits=[(w>>(7-i))&1 for i in range(8)]; return int(''.join(str(bits[p[i]]) for i in range(8)),2)
rng = np.random.default_rng(1)
for _ in range(200000):
    p = rng.permutation(8)
    Hp = [perm_word(w, p) for w in H]
    if set(H) & set(Hp) != {0, 255}: continue
    code = set()
    for x in Hp:
        for a in H:
            for b in H:
                code.add(((a ^ x) << 16) | ((b ^ x) << 8) | (a ^ b ^ x))
    if len(code) != 4096: continue
    wts = sorted(set(bin(c).count('1') for c in code))
    if wts == [0, 8, 12, 16, 24]:
        from collections import Counter
        print("perm", list(p)); print("H ", H); print("H'", Hp)
        print("weights", sorted(Counter(bin(c).count('1') for c in code).items()))
        np.save("turyn_H.npy", np.array([H, Hp]))
        break
