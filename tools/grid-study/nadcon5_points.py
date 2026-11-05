# NADCON5's source points against PROJ's grid: does the published grid come
# from NGS's thinned points (value at each point, converted to metres, against
# the point's), and what do the points cost to store?
# Usage: python3 -I nadcon5_points.py WORK      (after setup.sh, fetch_nadcon5.sh)
import json, subprocess, sys, numpy as np
W = sys.argv[1]; D = f'{W}/nadcon5-download'; T = 'nad83_2007.nad83_2011.conus'
def load(f): return np.array([[float(x) for x in l.split()[:3]] for l in open(f)])
pts = {b: load(f'{D}/smtcd{b}.{T}.60') for b in ('lat', 'lon', 'eht')}
for b, p in pts.items(): print(f'{b}: {len(p)} points, values {p[:,2].min():+.3f}..{p[:,2].max():+.3f} m')
print('lat and lon share positions:', np.array_equal(pts['lat'][:, :2], pts['lon'][:, :2]))
g_ = 'us_noaa_nadcon5_nad83_2007_nad83_2011_conus'
j = json.load(open(f'{W}/grids/{g_}.json')); w, h = j['size']; t = j['geoTransform']
g = np.fromfile(f'{W}/grids/{g_}.bin', dtype='<f4').reshape(3, h, w).astype(float)
for bi, b in enumerate(('lat', 'lon', 'eht')):
    p = pts[b]; lon = np.where(p[:, 0] > 180, p[:, 0] - 360, p[:, 0]); lat = p[:, 1]
    x = (lon - t[0]) / t[1] - 0.5; y = (lat - t[3]) / t[5] - 0.5
    x0 = np.floor(x).astype(int); y0 = np.floor(y).astype(int); fx = x - x0; fy = y - y0
    ok = (x0 >= 0) & (y0 >= 0) & (x0 < w - 1) & (y0 < h - 1)
    gv = (g[bi, y0, x0] * (1 - fx) * (1 - fy) + g[bi, y0, x0 + 1] * fx * (1 - fy)
          + g[bi, y0 + 1, x0] * (1 - fx) * fy + g[bi, y0 + 1, x0 + 1] * fx * fy)[ok]
    pv = p[ok, 2]; k = np.polyfit(gv, pv, 1)[0]
    print(f'  {b}: grid at the points vs point values: correlation {np.corrcoef(gv, pv)[0,1]:.4f}, '
          f'scale point/grid {k:.4g} (m per arc-second, or 1 for height), median |residual| {np.median(np.abs(pv - k*gv)):.4g} m')
def zstd(b): return len(subprocess.run(['zstd', '-19', '-q', '-c'], input=b, capture_output=True).stdout)
zz = lambda r: ((r << 1) ^ (r >> 63)).astype(np.uint32).view(np.uint8).reshape(-1, 4).T.copy().tobytes()
buf = b''
for b in ('lat', 'eht'):   # lon shares lat's positions
    p = pts[b]; o = np.lexsort((p[:, 0], p[:, 1])); pos = np.rint(p[o, :2] * 1e7).astype(np.int64)
    buf += zz(np.diff(pos, axis=0, prepend=0).ravel())
for b in ('lat', 'lon', 'eht'):
    ref = pts['eht' if b == 'eht' else 'lat']; o = np.lexsort((ref[:, 0], ref[:, 1]))
    buf += zz(np.rint(pts[b][o, 2] * 1000).astype(np.int64))
print(f'source points, all three bands (positions to 1e-7 deg, values to the mm as written): {zstd(buf)/1e6:.3f} MB')
