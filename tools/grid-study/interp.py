# Error-bounded encoding, SZ3-style: quantize to q, then predict each value
# coarse-to-fine by cubic interpolation of already-coded values, entropy-code
# the residuals. Lossless on the quantized field, so max error stays q/2.
import json, os, sys, subprocess, numpy as np
d = sys.argv[1]
def zstd(buf):
    p = os.path.join(d, 'tmp_i'); open(p, 'wb').write(buf)
    n = len(subprocess.run(['zstd', '-19', '--long=27', '-q', '-c', p], capture_output=True).stdout)
    os.remove(p); return n
def shuffle(a):
    return a.view(np.uint8).reshape(-1, a.itemsize).T.copy().tobytes()
def interp1(v, axis):
    # predict odd entries of v along axis from even ones: cubic where 4 neighbours exist, else linear/copy
    v = np.moveaxis(v, axis, 0)
    ev = v[0::2]; n_odd = v[1::2].shape[0]
    pred = np.empty((n_odd,) + v.shape[1:], dtype=np.float64)
    for k in range(n_odd):
        a, b = k, k + 1
        if b < ev.shape[0]:
            if a - 1 >= 0 and b + 1 < ev.shape[0]:
                pred[k] = (-ev[a - 1] + 9 * ev[a] + 9 * ev[b] - ev[b + 1]) / 16
            else:
                pred[k] = (ev[a] + ev[b]) / 2
        else:
            pred[k] = ev[a]
    return np.moveaxis(np.rint(pred).astype(np.int64), 0, axis)
def encode(iq):
    h, w = iq.shape; out = []
    L = int(np.floor(np.log2(max(h, w)))) - 3
    s = 2 ** max(L, 0)
    out.append(iq[::s, ::s].ravel())          # coarsest level verbatim (tiny)
    while s > 1:
        hs = s // 2
        sub = iq[::hs, ::hs]                  # this level's lattice; even/even known
        known = sub[0::2, 0::2]
        # rows: fill (even row, odd col)
        r = sub[0::2, :]
        p = interp1(r.astype(np.float64), 1)
        out.append((r[:, 1::2] - p).ravel())
        # cols: fill (odd row, every col) from even rows, now complete
        p = interp1(sub.astype(np.float64), 0)
        out.append((sub[1::2, :] - p).ravel())
        s = hs
    r = np.concatenate(out)
    return ((r << 1) ^ (r >> 63)).astype(np.uint32)
for g in sys.argv[2:]:
    j = json.load(open(os.path.join(d, g + '.json'))); w, h = j['size']; bands = j['bands']
    raw = np.fromfile(os.path.join(d, g + '.bin'), dtype='<f4').reshape(len(bands), h, w)
    for q_h, q_v in ((1e-5, 1e-4), (1e-4, 1e-3)):
        parts = []
        for b, band in enumerate(bands):
            desc = band.get('description', '')
            if 'accuracy' in desc: continue
            q = q_h if ('offset' in desc and 'height' not in desc) else q_v
            a = raw[b].astype(np.float64); mask = ~np.isfinite(a); a = np.where(mask, 0, a)
            iq = np.rint(a / q).astype(np.int64)
            parts.append(shuffle(encode(iq)))
            if mask.any(): parts.append(np.packbits(mask).tobytes())
        print(f"{g:45} q={q_h:g}/{q_v:g}: {zstd(b''.join(parts))/1e6:7.3f} MB")
