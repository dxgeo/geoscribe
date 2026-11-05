# EGM2008 as coefficients, sized with a degree-dependent step (the fix egm_errcheck.py
# found), over the band |lat| <= 57 deg of
# PROJ's 2.5' grid (the sectoral recursion is scaled by 1e280 and computed in
# logs, which holds to degree 2190 below about 57 deg).
import json, os, sys, subprocess, time, numpy as np
S = sys.argv[1]; LATMAX = 57.0
z = np.load(f'{S}/egm-work/egm2008.npz'); C, Sn, N = z['C'], z['S'], int(z['N'])
GM, a = 0.3986004415e15, 0.63781363e07
GMw, aw, e2 = 3.986004418e14, 6378137.0, 0.00669437999013
for n, c in ((2, -0.484166774985e-3), (4, 0.790303733511e-6), (6, -0.168724961151e-8), (8, 0.346052468394e-11), (10, -0.265002225747e-14)):
    C[n, 0] -= c * (GMw / GM) * (aw / a) ** n
C[0, 0] = 0.0
def zstd(buf):
    p = f'{S}/egm-work/tmp8'; open(p, 'wb').write(buf)
    k = len(subprocess.run(['zstd', '-19', '--long=27', '-T0', '-q', '-c', p], capture_output=True).stdout); os.remove(p); return k
def zz(r):
    r = r.astype(np.int64); return ((r << 1) ^ (r >> 63)).astype(np.uint32)
def shuffle(a): return a.view(np.uint8).reshape(-1, a.itemsize).T.copy().tobytes()
j = json.load(open(f'{S}/grids/us_nga_egm08_25.json')); w, h = j['size']; t = j['geoTransform']
full = np.fromfile(f'{S}/grids/us_nga_egm08_25.bin', dtype='<f4').reshape(h, w).astype(np.float64)
lat_all = t[3] + (np.arange(h) + 0.5) * t[5]; lon = t[0] + (np.arange(w) + 0.5) * t[1]
rows = np.where(np.abs(lat_all) <= LATMAX)[0]; grid = full[rows]; lat = lat_all[rows]; hb = len(rows)
phi = np.radians(lat); sp = np.sin(phi)
Nv = aw / np.sqrt(1 - e2 * sp ** 2); x = Nv * np.cos(phi); zc = Nv * (1 - e2) * sp
r = np.hypot(x, zc); tt = zc / r; uu = x / r
gamma = 9.7803253359 * (1 + 0.00193185265241 * sp ** 2) / np.sqrt(1 - e2 * sp ** 2)
q = a / r; lnq = np.log(q); lnu = np.log(uu)
SCALE = 1e280
def synth(C, Sn):
    A = np.zeros((hb, N + 1)); B = np.zeros((hb, N + 1))
    lnpmm = np.zeros(hb)
    for m in range(N + 1):
        if m == 1: lnpmm = np.log(np.sqrt(3.0)) + lnu
        elif m > 1: lnpmm = lnpmm + lnu + 0.5 * np.log((2 * m + 1) / (2 * m))
        pmm = np.exp(np.maximum(lnpmm + np.log(SCALE), -745.0))   # scaled
        pmm[lnpmm + np.log(SCALE) < -740] = 0.0
        p2 = np.zeros(hb); p1 = pmm
        cm = C[:, m]; sm = Sn[:, m]
        for n in range(m, N + 1):
            if n == m: p = pmm
            elif n == m + 1: p = np.sqrt(2 * m + 3.0) * tt * pmm
            else:
                anm = np.sqrt((2 * n - 1.0) * (2 * n + 1) / ((n - m) * (n + m)))
                bnm = np.sqrt((2 * n + 1.0) * (n + m - 1) * (n - m - 1) / ((n - m) * (n + m) * (2 * n - 3)))
                p = anm * tt * p1 - bnm * p2
            if n >= m + 1: p2, p1 = p1, p
            if n >= 2 and (cm[n] != 0.0 or sm[n] != 0.0):
                f = np.exp(n * lnq) * p
                A[:, m] += f * cm[n]; B[:, m] += f * sm[n]
    A /= SCALE; B /= SCALE
    lam = np.radians(lon); ms = np.arange(N + 1)
    T = A @ np.cos(np.outer(ms, lam)) + B @ np.sin(np.outer(ms, lam))
    return (GM / (r * gamma))[:, None] * T

# degree-dependent step: d_n = d0 * (r_edge / a)^n, so the (a/r)^n growth at the
# band's edge (or the pole) is paid for in the step rather than in the error
tri = np.tril_indices(N + 1); ncoef = (N + 1) * (N + 2) - (N + 1); keep = tri[1] > 0
def r_at(lat_deg):
    p_ = np.radians(lat_deg); s_ = np.sin(p_); Nv_ = aw / np.sqrt(1 - e2 * s_ ** 2)
    return np.hypot(Nv_ * np.cos(p_), Nv_ * (1 - e2) * s_)
nn = np.arange(N + 1)[:, None] * np.ones((1, N + 1))
for edge, label in ((LATMAX, 'band edge 57 deg'), (90.0, 'pole (global grid)')):
    grow = (a / r_at(edge)) ** nn
    for target in (1e-3,):
        d0 = target * np.sqrt(12) / (a * np.sqrt(ncoef))
        d = d0 / grow
        iC = np.rint(C[tri] / d[tri]); iS = np.rint(Sn[tri] / d[tri])
        size = zstd(shuffle(zz(iC)) + shuffle(zz(iS[keep])))
        msg = f'  step for {target*1000:g} mm rms at the {label}: {size/1e6:.3f} MB'
        if edge == LATMAX:
            dC = np.zeros_like(C); dS = np.zeros_like(Sn); dC[tri] = iC * d[tri] - C[tri]; dS[tri] = iS * d[tri] - Sn[tri]
            err = synth(dC, dS)
            msg += f'; geoid change over the band rms {np.sqrt(np.mean(err**2))*1000:.3f} mm, max {np.max(np.abs(err))*1000:.3f} mm'
        print(msg, flush=True)
