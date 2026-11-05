# How far is zstd from the best an entropy coder could do on the same
# residuals? Order-0 entropy, and entropy conditioned on a context (the
# magnitude bucket of the left and upper residuals), as a JPEG-LS/CALIC-style
# coder would use. Bytes = bits * count / 8.
import json, os, sys, subprocess, numpy as np
d = sys.argv[1]; q_h = float(sys.argv[2]); q_v = float(sys.argv[3])
def H(x):
    _, c = np.unique(x, return_counts=True); p = c / c.sum(); return -(p * np.log2(p)).sum()
def ctxH(r):
    # context: bucket of |left|+|up| residual
    m = np.abs(np.pad(r, ((1, 0), (1, 0))))
    c = np.minimum(np.floor(np.log2(1 + m[1:, :-1] + m[:-1, 1:])).astype(int), 15)
    tot = 0.0
    for k in np.unique(c):
        sel = r[c == k]; tot += H(sel) * sel.size
    return tot / r.size
def zstd(buf):
    p = os.path.join(d, 'tmp_e'); open(p, 'wb').write(buf)
    n = len(subprocess.run(['zstd', '-19', '--long=27', '-q', '-c', p], capture_output=True).stdout); os.remove(p); return n
for g in sys.argv[4:]:
    j = json.load(open(os.path.join(d, g + '.json'))); w, h = j['size']; bands = j['bands']
    raw = np.fromfile(os.path.join(d, g + '.bin'), dtype='<f4').reshape(len(bands), h, w)
    z = h0 = hc = 0
    for b, band in enumerate(bands):
        desc = band.get('description', '')
        if 'accuracy' in desc or desc == '' and b >= 2: continue
        q = q_h if ('offset' in desc and 'height' not in desc) else q_v
        a = np.nan_to_num(raw[b].astype(np.float64)); iq = np.rint(a / q).astype(np.int64); r = iq.copy()
        r[1:, 1:] = iq[1:, 1:] - iq[:-1, 1:] - iq[1:, :-1] + iq[:-1, :-1]; r[0, 1:] = iq[0, 1:] - iq[0, :-1]; r[1:, 0] = iq[1:, 0] - iq[:-1, 0]
        zz = ((r << 1) ^ (r >> 63)).astype(np.uint32)
        z += zstd(zz.view(np.uint8).reshape(-1, 4).T.copy().tobytes())
        h0 += H(r) * r.size / 8; hc += ctxH(r) * r.size / 8
    print(f"{g:44} zstd {z/1e6:7.3f} MB | order-0 entropy {h0/1e6:7.3f} MB | context entropy {hc/1e6:7.3f} MB")
