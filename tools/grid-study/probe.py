# For each grid band: (1) its native decimal step -- the coarsest 10^-k step
# that reproduces every value to float32 rounding; (2) roughness: the RMS of
# the 2D second difference, in the band's unit, against a pure-noise model;
# (3) upsampling: whether interior nodes are an exact interpolation of a
# coarser lattice (fraction of near-zero second differences per stride).
import json, os, sys, numpy as np
d = sys.argv[1]
for g in sys.argv[2:]:
    j = json.load(open(os.path.join(d, g + '.json'))); w, h = j['size']; bands = j['bands']
    raw = np.fromfile(os.path.join(d, g + '.bin'), dtype='<f4').reshape(len(bands), h, w)
    for b, band in enumerate(bands):
        desc = band.get('description', ''); unit = band.get('unit', '')
        if 'accuracy' in desc: continue
        a = raw[b].astype(np.float64); a = a[np.isfinite(a)].reshape(-1) if not np.isfinite(raw[b]).all() else raw[b].astype(np.float64)
        flat = a.ravel()
        native = None
        for k in range(1, 9):
            q = 10.0 ** -k
            err = np.max(np.abs(np.rint(flat / q) * q - flat))
            tol = np.max(np.abs(flat)) * 2 ** -23 * 2
            if err <= max(tol, 1e-12): native = q; break
        line = f"{g[:34]:34} {desc[:18]:18} native step: {native if native else 'none (full float32)'}"
        if a.ndim == 2:
            lap = a[1:-1, 1:-1] * 4 - a[:-2, 1:-1] - a[2:, 1:-1] - a[1:-1, :-2] - a[1:-1, 2:]
            rms = np.sqrt(np.mean(lap ** 2))
            scale = 30.9 if 'arc-second' in unit else 1.0   # arcsec -> metres (latitude)
            # upsampling: share of rows whose second difference along x is ~0, per stride
            d2 = a[:, 2:] - 2 * a[:, 1:-1] + a[:, :-2]
            col_rms = np.sqrt(np.mean(d2 ** 2, axis=0))
            per = {s: float(np.mean(col_rms[np.arange(len(col_rms)) % s != s - 1] < 1e-3 * np.median(col_rms[np.arange(len(col_rms)) % s == s - 1] + 1e-30))) for s in (2, 3, 4, 5, 8, 10)}
            best = max(per, key=per.get)
            line += f" | roughness {rms*scale*1000:8.3f} mm | upsampled-by-{best}? {per[best]:.0%} of in-between columns linear"
        print(line)
