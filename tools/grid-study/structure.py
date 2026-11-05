# (1) a lattice in metres: latitude offsets (arcsec) converted to metres along
#     the meridian at each row; how close to whole multiples of 1e-3 / 1e-4 m?
# (2) filler: share of nodes that are no-data, exactly zero, or in long runs of
#     one repeated value, and how the residual cost splits by the field's activity.
import json, os, sys, numpy as np
d = sys.argv[1]; a, e2 = 6378137.0, 0.00669437999013
for g in sys.argv[2:]:
    j = json.load(open(os.path.join(d, g + '.json'))); w, h = j['size']; t = j['geoTransform']; bands = j['bands']
    raw = np.fromfile(os.path.join(d, g + '.bin'), dtype='<f4').reshape(len(bands), h, w).astype(np.float64)
    lat = t[3] + (np.arange(h) + 0.5) * t[5]
    print(f'{g}  {w}x{h}  spacing {abs(t[1])*3600:.3g}" x {abs(t[5])*3600:.3g}"  lat {lat.min():.2f}..{lat.max():.2f}')
    for b, band in enumerate(bands[:2]):
        desc = band.get('description', ''); v = raw[b]
        nd = band.get('noDataValue'); bad = ~np.isfinite(v) | ((v == nd) if nd is not None else False)
        zero = (v == 0) & ~bad
        # runs of an identical value along rows
        same = np.zeros_like(v, dtype=bool); same[:, 1:] = v[:, 1:] == v[:, :-1]
        line = f'   {desc[:17]:17} no-data {bad.mean():6.1%}  exact zero {zero.mean():6.1%}  repeats left neighbour {same.mean():6.1%}'
        if desc == 'latitude_offset' and band.get('unit') == 'arc-second':
            phi = np.radians(lat); M = a * (1 - e2) / (1 - e2 * np.sin(phi) ** 2) ** 1.5
            metres = v * np.pi / 648000 * M[:, None]
            for step in (1e-3, 1e-4):
                frac = np.abs(metres / step - np.rint(metres / step))
                line += f' | metres on a {step*1000:g} mm lattice: {np.mean(frac[~bad] < 0.01):.1%} (chance 2%)'
        print(line)
    # where does the cost go: Lorenzo residual magnitude vs local activity
    v = np.nan_to_num(raw[0]); iq = np.rint(v / 1e-5)
    r = iq[1:, 1:] - iq[:-1, 1:] - iq[1:, :-1] + iq[:-1, :-1]
    bits = np.log2(1 + np.abs(r)) + 1
    flat = np.abs(r) <= 1
    print(f'   residual (0.3 mm step): {flat.mean():.1%} of nodes cost about nothing (|r|<=1); they hold {bits[flat].sum()/bits.sum():.1%} of the bits')
