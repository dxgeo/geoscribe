# PROJ's NADCON5 grid against GMT 6.7 surface run with NGS's published command
# on NGS's published thinned points, node by node, in millimetres.
import json, sys, numpy as np
S = sys.argv[1]; W = f'{S}/gmt-work'   # Usage: python3 -I cmp_surface.py WORK
j = json.load(open(f'{S}/grids/us_noaa_nadcon5_nad83_2007_nad83_2011_conus.json')); w, h = j['size']; t = j['geoTransform']
proj = np.fromfile(f'{S}/grids/us_noaa_nadcon5_nad83_2007_nad83_2011_conus.bin', dtype='<f4').reshape(3, h, w).astype(np.float64)
lat = t[3] + (np.arange(h) + 0.5) * t[5]
a, e2 = 6378137.0, 0.00669438002290          # GRS80
s = np.sin(np.radians(lat))[:, None]
M = a * (1 - e2) / (1 - e2 * s ** 2) ** 1.5; Nn = a / np.sqrt(1 - e2 * s ** 2)
to_m = [np.pi / 648000 * M, np.pi / 648000 * Nn * np.cos(np.radians(lat))[:, None], 1.0]
for i, b in enumerate(('lat', 'lon', 'eht')):
    g = np.fromfile(f'{W}/{b}.bin', dtype='<f8').reshape(h, w)
    p = proj[i] * to_m[i]
    d = (g - p) * 1000
    pc = np.percentile(np.abs(d), [50, 90, 99, 99.9, 100])
    print(f'{b}: PROJ minus GMT 6.7, |mm|: p50 {pc[0]:.3f}  p90 {pc[1]:.3f}  p99 {pc[2]:.3f}  p99.9 {pc[3]:.3f}  max {pc[4]:.3f}   (field spread {np.std(p)*1000:.1f} mm)')
    if b == 'lat':
        rows = np.array_split(np.arange(h), 4); cols = np.array_split(np.arange(w), 6)
        print('   median |diff| by region (mm), north at top:')
        for r in rows: print('   ' + ' '.join(f'{np.median(np.abs(d[np.ix_(r, c)])):7.3f}' for c in cols))

