# EGM96 as coefficients: (1) what the coefficients cost, quantized so the geoid
# they give moves by a bounded amount -- checked by synthesis, not estimated;
# (2) how much of PROJ's 15' EGM96 grid the model explains, and what the rest
# costs to store.
import json, os, sys, subprocess, numpy as np
S = sys.argv[1]
z = np.load(f'{S}/egm-work/egm96.npz'); C, Sn, sC, sS, N = z['C'], z['S'], z['sC'], z['sS'], int(z['N'])
GM, a = 0.3986004415e15, 0.6378136300e07
# WGS84 normal field even zonals (fully normalised), rescaled to EGM96's GM and a
GMw, aw, e2 = 3.986004418e14, 6378137.0, 0.00669437999013
for n, c in ((2, -0.484166774985e-3), (4, 0.790303733511e-6), (6, -0.168724961151e-8), (8, 0.346052468394e-11), (10, -0.265002225747e-14)):
    C[n, 0] -= c * (GMw / GM) * (aw / a) ** n
C[0, 0] = 0.0
def zstd(buf):
    p = f'{S}/egm-work/tmp'; open(p, 'wb').write(buf)
    k = len(subprocess.run(['zstd', '-19', '--long=27', '-q', '-c', p], capture_output=True).stdout); os.remove(p); return k
def zz(r):
    r = r.astype(np.int64); return ((r << 1) ^ (r >> 63)).astype(np.uint32)
def shuffle(a): return a.view(np.uint8).reshape(-1, a.itemsize).T.copy().tobytes()
# the grid PROJ ships
j = json.load(open(f'{S}/grids/us_nga_egm96_15.json')); w, h = j['size']; t = j['geoTransform']
grid = np.fromfile(f'{S}/grids/us_nga_egm96_15.bin', dtype='<f4').reshape(h, w).astype(np.float64)
lat = t[3] + (np.arange(h) + 0.5) * t[5]; lon = t[0] + (np.arange(w) + 0.5) * t[1]
print('grid', w, 'x', h, 'lat', lat[0], '..', lat[-1], 'lon', lon[0], '..', lon[-1])
phi = np.radians(lat); sp = np.sin(phi)
Nv = aw / np.sqrt(1 - e2 * sp ** 2); x = Nv * np.cos(phi); zc = Nv * (1 - e2) * sp
r = np.hypot(x, zc); tt = zc / r; uu = x / r
gamma = 9.7803253359 * (1 + 0.00193185265241 * sp ** 2) / np.sqrt(1 - e2 * sp ** 2)
q = a / r
def synth(C, Sn):
    A = np.zeros((h, N + 1)); B = np.zeros((h, N + 1))
    pmm = np.ones(h)
    for m in range(N + 1):
        if m == 1: pmm = np.sqrt(3.0) * uu
        elif m > 1: pmm = pmm * uu * np.sqrt((2 * m + 1) / (2 * m))
        p2 = np.zeros(h); p1 = pmm
        for n in range(m, N + 1):
            if n == m: p = pmm
            elif n == m + 1: p = np.sqrt(2 * m + 3.0) * tt * pmm
            else:
                anm = np.sqrt((2 * n - 1.0) * (2 * n + 1) / ((n - m) * (n + m)))
                bnm = np.sqrt((2 * n + 1.0) * (n + m - 1) * (n - m - 1) / ((n - m) * (n + m) * (2 * n - 3)))
                p = anm * tt * p1 - bnm * p2
            if n >= m + 1: p2, p1 = p1, p
            if n >= 2:
                f = q ** n * p
                A[:, m] += f * C[n, m]; B[:, m] += f * Sn[n, m]
    lam = np.radians(lon); ms = np.arange(N + 1)
    T = A @ np.cos(np.outer(ms, lam)) + B @ np.sin(np.outer(ms, lam))
    return (GM / (r * gamma))[:, None] * T
zeta = synth(C, Sn)
res = grid - zeta
pc = np.percentile(np.abs(res - np.median(res)), [50, 90, 99, 100])
print(f'grid - model: median {np.median(res):+.3f} m; |deviation| p50 {pc[0]:.3f}, p90 {pc[1]:.3f}, p99 {pc[2]:.3f}, max {pc[3]:.3f} m')
def enc_grid(a, qs):
    iq = np.rint(a / qs).astype(np.int64); rr = iq.copy()
    rr[1:, 1:] = iq[1:, 1:] - iq[:-1, 1:] - iq[1:, :-1] + iq[:-1, :-1]; rr[0, 1:] = iq[0, 1:] - iq[0, :-1]; rr[1:, 0] = iq[1:, 0] - iq[:-1, 0]
    return zstd(shuffle(zz(rr)))
for qs in (1e-5, 1e-3):
    print(f'  step {qs*1000:g} mm: grid alone {enc_grid(grid, qs)/1e6:.3f} MB, residual from model {enc_grid(res, qs)/1e6:.3f} MB')
# the coefficients
tri = np.tril_indices(N + 1)
def coef_cost(dC, dS, label):
    iC = np.rint(C[tri] / dC[tri] if np.ndim(dC) else C[tri] / dC); iS = np.rint(Sn[tri] / dS[tri] if np.ndim(dS) else Sn[tri] / dS)
    keep = tri[1] > 0
    payload = shuffle(zz(iC)) + shuffle(zz(iS[keep]))
    Cq = np.zeros_like(C); Sq = np.zeros_like(Sn)
    Cq[tri] = iC * (dC[tri] if np.ndim(dC) else dC); Sq[tri] = iS * (dS[tri] if np.ndim(dS) else dS)
    err = synth(Cq, Sq) - zeta
    print(f'  {label:44} {zstd(payload)/1e6:.3f} MB; geoid change rms {np.sqrt(np.mean(err**2))*1000:.3f} mm, max {np.max(np.abs(err))*1000:.3f} mm')
ncoef = (N + 1) * (N + 2) - (N + 1)
print('coefficients:', ncoef, 'numbers;', f'{ncoef*8/1e6:.2f} MB as float64')
for target in (1e-3, 1e-4):
    d = target * np.sqrt(12) / (a * np.sqrt(ncoef))
    coef_cost(d, d, f'uniform step for {target*1000:g} mm rms')
sig = lambda s: np.where(s > 0, s, np.min(s[s > 0]))
coef_cost(sig(sC) / 4, sig(sS) / 4, 'step sigma/4 (a quarter of each one\'s error)')
