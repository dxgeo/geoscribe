# Sample grids, quantised to a step and neighbour-predicted, then zstd -19 and
# xz -9e: the first per-file measurements (plans/transformation-data.org
# § MEASURED). Steps by env: QH (arc-second offsets, default 1e-5 = 0.3 mm),
# QV (metres, default 1e-4 = 0.1 mm); DROP_ACC=1 drops accuracy bands.
# Usage: PROJ_DATA=... python3 -I gridpack.py WORK/grids GRID [GRID ...]
import json, subprocess, sys, os, numpy as np
d = sys.argv[1]
def comp(buf, tag):
    p = os.path.join(d, 'tmp_' + tag)
    open(p, 'wb').write(buf)
    z = subprocess.run(['zstd', '-19', '--long=27', '-q', '-c', p], capture_output=True).stdout
    x = subprocess.run(['xz', '-9e', '-c', p], capture_output=True).stdout
    os.remove(p)
    return len(z), len(x)
def shuffle(a):  # byte planes
    return a.view(np.uint8).reshape(-1, a.itemsize).T.copy().tobytes()
P = os.environ['PROJ_DATA']
print(f"{'grid':45} {'tif':>8} {'q+pred zstd':>11} {'q+pred xz':>10}  quantum / max err")
for g in sys.argv[2:]:
    j = json.load(open(os.path.join(d, g + '.json')))
    w, h = j['size']; bands = j['bands']
    raw = np.fromfile(os.path.join(d, g + '.bin'), dtype='<f4').reshape(len(bands), h, w)
    tif = os.path.getsize(os.path.join(P, g + '.tif')) if os.path.exists(os.path.join(P, g + '.tif')) else 0
    parts = []; errs = []
    for b, band in enumerate(bands):
        desc = band.get('description', '')
        q = float(os.environ.get('QH','1e-5')) if 'offset' in desc and 'accuracy' not in desc and 'height' not in desc else float(os.environ.get('QV','1e-4'))
        if 'accuracy' in desc:
            if os.environ.get('DROP_ACC'): continue
            q = 1e-3
        a = raw[b].astype(np.float64)
        mask = ~np.isfinite(a)
        a = np.where(mask, 0.0, a)
        iq = np.round(a / q).astype(np.int64)
        errs.append((desc[:20], q, np.max(np.abs(iq * q - a))))
        r = iq.copy()
        r[1:, 1:] = iq[1:, 1:] - iq[:-1, 1:] - iq[1:, :-1] + iq[:-1, :-1]
        r[0, 1:] = iq[0, 1:] - iq[0, :-1]
        r[1:, 0] = iq[1:, 0] - iq[:-1, 0]
        zz = ((r << 1) ^ (r >> 63)).astype(np.uint32)
        parts.append(shuffle(zz.ravel()))
        if mask.any(): parts.append(np.packbits(mask).tobytes())
    z, x = comp(b''.join(parts), 'q')
    e = '; '.join(f"{n}:{q:g}/{m:.1e}" for n, q, m in errs[:2])
    print(f"{g:45} {tif/1e6:7.2f}M {z/1e6:10.3f}M {x/1e6:9.3f}M  {e}")
