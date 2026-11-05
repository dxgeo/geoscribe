# Cross-grid prediction: encode a national geoid as its residual from the
# global model (EGM2008, bilinear at the national grid's nodes), against the
# national geoid encoded alone. Same quantize + neighbour predictor + zstd.
import json, os, sys, subprocess, numpy as np
d = sys.argv[1]
def load(g):
    j = json.load(open(os.path.join(d, g + '.json'))); w, h = j['size']
    a = np.fromfile(os.path.join(d, g + '.bin'), dtype='<f4').reshape(len(j['bands']), h, w)[0].astype(np.float64)
    nd = j['bands'][0].get('noDataValue')
    if nd is not None: a[a == nd] = np.nan
    return a, j['geoTransform']
def zstd(buf):
    p = os.path.join(d, 'tmp_x'); open(p, 'wb').write(buf)
    n = len(subprocess.run(['zstd', '-19', '--long=27', '-q', '-c', p], capture_output=True).stdout); os.remove(p); return n
def enc(a, q):
    mask = ~np.isfinite(a); a = np.where(mask, 0, a)
    iq = np.rint(a / q).astype(np.int64); r = iq.copy()
    r[1:, 1:] = iq[1:, 1:] - iq[:-1, 1:] - iq[1:, :-1] + iq[:-1, :-1]
    r[0, 1:] = iq[0, 1:] - iq[0, :-1]; r[1:, 0] = iq[1:, 0] - iq[:-1, 0]
    zz = ((r << 1) ^ (r >> 63)).astype(np.uint32)
    buf = zz.view(np.uint8).reshape(-1, 4).T.copy().tobytes() + (np.packbits(mask).tobytes() if mask.any() else b'')
    return zstd(buf) / 1e6
E, et = load('us_nga_egm08_25')
def egm_at(lon, lat):
    x = (lon - et[0]) / et[1] - 0.5; y = (lat - et[3]) / et[5] - 0.5   # pixel-is-area centres
    x = np.mod(x, E.shape[1]); x0 = np.floor(x).astype(int); y0 = np.clip(np.floor(y).astype(int), 0, E.shape[0] - 2)
    fx = x - x0; fy = y - y0; x1 = (x0 + 1) % E.shape[1]
    return (E[y0, x0] * (1 - fx) * (1 - fy) + E[y0, x1] * fx * (1 - fy) + E[y0 + 1, x0] * (1 - fx) * fy + E[y0 + 1, x1] * fx * fy)
for g in sys.argv[2:]:
    a, t = load(g); h, w = a.shape
    off = 0.5
    lon = t[0] + (np.arange(w) + off) * t[1]; lat = t[3] + (np.arange(h) + off) * t[5]
    LON, LAT = np.meshgrid(lon, lat)
    res = a - egm_at(LON, LAT)
    ok = np.isfinite(res)
    tif = os.path.getsize(os.path.join(os.environ['PROJ_DATA'], g + '.tif')) / 1e6
    print(f"{g:30} {w}x{h} tif {tif:6.2f} MB | residual from EGM2008: rms {np.sqrt(np.nanmean(res**2)):.3f} m, range {np.nanmin(res):+.2f}..{np.nanmax(res):+.2f} m")
    for q in (1e-5, 1e-4, 1e-3):
        print(f"    q={q*1000:g} mm: alone {enc(a, q):7.3f} MB   as residual {enc(res, q):7.3f} MB")
