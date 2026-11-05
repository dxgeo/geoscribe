# Parse an ICGEM .gfc file into dense lower-triangular C, S, sigma arrays (.npz).
import sys, numpy as np
src, dst = sys.argv[1], sys.argv[2]
with open(src) as f:
    for line in f:
        if line.startswith('end_of_head'): break
        if line.startswith('max_degree'): N = int(line.split()[1])
    body = f.read().replace('d', 'e').replace('D', 'E')
tok = np.array(body.split(), dtype=object).reshape(-1, 7)
assert (tok[:, 0] == 'gfc').all()
n = tok[:, 1].astype(int); m = tok[:, 2].astype(int)
C = np.zeros((N + 1, N + 1)); S = np.zeros((N + 1, N + 1)); sC = np.zeros((N + 1, N + 1)); sS = np.zeros((N + 1, N + 1))
C[n, m] = tok[:, 3].astype(float); S[n, m] = tok[:, 4].astype(float); sC[n, m] = tok[:, 5].astype(float); sS[n, m] = tok[:, 6].astype(float)
np.savez(dst, C=C, S=S, sC=sC, sS=sS, N=N)
print(src.split('/')[-1], 'max degree', N, 'coefficient lines', len(n))
