# GMT 6.7's rebuild of NADCON5 (nadcon5_surface.sh) against PROJ's grid, by
# each node's distance to the nearest thinned control point -- in grid cells
# (1' ~ 1.85 km north-south, ~1.4 km east-west), grown out from the cells that
# hold a point, a cell a step.
# Usage: python3 -I nadcon5_distance.py WORK
import json, sys, numpy as np
W = sys.argv[1]; D = f'{W}/nadcon5-download'; G = f'{W}/gmt-work'; T = 'nad83_2007.nad83_2011.conus'
g_ = 'us_noaa_nadcon5_nad83_2007_nad83_2011_conus'
j = json.load(open(f'{W}/grids/{g_}.json')); w, h = j['size']; t = j['geoTransform']
proj = np.fromfile(f'{W}/grids/{g_}.bin', dtype='<f4').reshape(3, h, w).astype(float)
lat = t[3] + (np.arange(h) + 0.5) * t[5]; a, e2 = 6378137.0, 0.00669438002290
s = np.sin(np.radians(lat))[:, None]; M = a * (1 - e2) / (1 - e2 * s ** 2) ** 1.5; Nn = a / np.sqrt(1 - e2 * s ** 2)
to_m = [np.pi / 648000 * M, np.pi / 648000 * Nn * np.cos(np.radians(lat))[:, None], 1.0]
p = np.array([[float(x) for x in l.split()[:2]] for l in open(f'{D}/smtcdlat.{T}.60')])
lon = np.where(p[:, 0] > 180, p[:, 0] - 360, p[:, 0])
c = np.rint((lon - t[0]) / t[1] - 0.5).astype(int); r = np.rint((p[:, 1] - t[3]) / t[5] - 0.5).astype(int)
ok = (c >= 0) & (c < w) & (r >= 0) & (r < h)
dist = np.full((h, w), 999, dtype=np.int32); cur = np.zeros((h, w), bool); cur[r[ok], c[ok]] = True; dist[cur] = 0
for k in range(1, 151):
    nxt = cur.copy()
    nxt[1:] |= cur[:-1]; nxt[:-1] |= cur[1:]; nxt[:, 1:] |= cur[:, :-1]; nxt[:, :-1] |= cur[:, 1:]
    nxt[1:, 1:] |= cur[:-1, :-1]; nxt[:-1, :-1] |= cur[1:, 1:]; nxt[1:, :-1] |= cur[:-1, 1:]; nxt[:-1, 1:] |= cur[1:, :-1]
    dist[nxt & ~cur] = k; cur = nxt
print('nodes by distance (cells) to the nearest control point:', ', '.join(f'<={k}: {np.mean(dist<=k):.0%}' for k in (1, 3, 6, 15, 30, 60)))
for i, b in enumerate(('lat', 'lon', 'eht')):
    g = np.fromfile(f'{G}/{b}.bin', dtype='<f8').reshape(h, w); d = np.abs(g - proj[i] * to_m[i]) * 1000
    print(f' {b}:')
    for lo, hi in ((0, 1), (2, 3), (4, 6), (7, 15), (16, 30), (31, 60), (61, 999)):
        sel = (dist >= lo) & (dist <= hi)
        if sel.any():
            print(f'   {lo:>3}-{hi:<3} cells ({sel.mean():5.1%} of nodes): |diff| p50 {np.median(d[sel]):7.3f}  p99 {np.percentile(d[sel],99):7.3f}  max {d[sel].max():7.3f} mm')
