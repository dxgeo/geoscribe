# Every grid file proj.db references and PROJ has installed: each subgrid,
# each band read band-sequential, accuracy bands dropped, quantised by the
# band's unit (arc-second, degree, metre), neighbour-predicted, zstd'd.
# Prints one line per file and the totals, at two steps.
import json, os, subprocess, sys, tempfile, numpy as np
S = sys.argv[1]; P = os.environ['PROJ_DATA']
STEPS = {'fine': {'arc-second': 1e-5, 'metre': 1e-4}, 'coarse': {'arc-second': 1e-4, 'metre': 1e-3}}
files = {}
for l in open(f'{S}/grids.txt'):
    m, g = l.strip().split('|', 1)
    if os.path.exists(os.path.join(P, g)) and g.endswith('.tif'):
        files[g] = 'vertical' if any(k in m.lower() for k in ('height', 'geoid', 'vertical', 'depth')) else 'horizontal'
def zstd_len(buf):
    r = subprocess.run(['zstd', '-19', '-T0', '--long=27', '-q', '-c'], input=buf, capture_output=True); return len(r.stdout)
def residual(iq):
    r = iq.copy()
    r[1:, 1:] = iq[1:, 1:] - iq[:-1, 1:] - iq[1:, :-1] + iq[:-1, :-1]; r[0, 1:] = iq[0, 1:] - iq[0, :-1]; r[1:, 0] = iq[1:, 0] - iq[:-1, 0]
    z = ((r << 1) ^ (r >> 63)).astype(np.uint32); return z.view(np.uint8).reshape(-1, 4).T.copy().tobytes()
tmp = tempfile.mkdtemp(dir=S)
tot = {}; skipped = []
done = {}  # resume: rows already measured, from a previous run's output
if os.path.exists(f'{S}/inventory.done'):
    for l in open(f'{S}/inventory.done'):
        k, g, sh, _, f, _, c, _ = l.split()
        done[g] = (k, float(sh) * 1e6, float(f) * 1e6, float(c) * 1e6)
for g, kind in sorted(files.items()):
    if g in done:
        k, sh, f, c = done[g]; t = tot.setdefault(k, [0, 0, 0, 0]); t[0] += 1; t[1] += sh; t[2] += f; t[3] += c
        print(f'{k:10} {g:60} {sh/1e6:8.2f} -> {f/1e6:7.3f} / {c/1e6:7.3f} MB', flush=True); continue
    info = json.loads(subprocess.run(['gdalinfo', '-json', os.path.join(P, g)], capture_output=True, text=True).stdout)
    subs = [v for k, v in info.get('metadata', {}).get('SUBDATASETS', {}).items() if k.endswith('_NAME')] or [os.path.join(P, g)]
    parts = {'fine': [], 'coarse': []}
    for sd in subs:
        j = json.loads(subprocess.run(['gdalinfo', '-json', sd], capture_output=True, text=True).stdout)
        w, h = j['size']; bands = j['bands']
        b = os.path.join(tmp, 'g.bin')
        subprocess.run(['gdal_translate', '-q', '-of', 'ENVI', '-co', 'INTERLEAVE=BSQ', '-ot', 'Float32', sd, b], check=True)
        raw = np.fromfile(b, dtype='<f4').reshape(len(bands), h, w)
        for i, band in enumerate(bands):
            desc = band.get('description', '') or ''; unit = band.get('unit', '') or ''
            if 'accuracy' in desc or (desc == '' and len(bands) == 4 and i >= 2): continue
            if unit == 'arc-second' or (unit == '' and 'offset' in desc and 'height' not in desc and kind == 'horizontal'): u = 'arc-second'
            elif unit == 'degree': u = 'degree'
            elif unit in ('metre', 'm', 'meter') or kind == 'vertical' or 'height' in desc: u = 'metre'
            else: skipped.append(f'{g}:{desc}:{unit}'); continue
            a = raw[i].astype(np.float64); nd = band.get('noDataValue')
            mask = ~np.isfinite(a) | ((a == nd) if nd is not None else False); a = np.where(mask, 0.0, a)
            for name, st in STEPS.items():
                q = st['arc-second'] / 3600 if u == 'degree' else st[u]
                parts[name].append(residual(np.rint(a / q).astype(np.int64)))
                if mask.any(): parts[name].append(np.packbits(mask).tobytes())
        os.remove(b)
    shipped = os.path.getsize(os.path.join(P, g))
    sizes = {name: zstd_len(b''.join(p)) for name, p in parts.items()}
    t = tot.setdefault(kind, [0, 0, 0, 0]); t[0] += 1; t[1] += shipped; t[2] += sizes['fine']; t[3] += sizes['coarse']
    print(f'{kind:10} {g:60} {shipped/1e6:8.2f} -> {sizes["fine"]/1e6:7.3f} / {sizes["coarse"]/1e6:7.3f} MB', flush=True)
print('\nTOTALS (shipped -> 0.3 mm & 0.1 mm step / 3 mm & 1 mm step)')
for kind, (n, s, f, c) in tot.items():
    print(f'  {kind:10} {n:3d} files  {s/1e6:8.1f} MB -> {f/1e6:7.1f} MB / {c/1e6:7.1f} MB   ({s/f:.1f}x / {s/c:.1f}x)')
print('skipped bands:', len(skipped), skipped[:10])
