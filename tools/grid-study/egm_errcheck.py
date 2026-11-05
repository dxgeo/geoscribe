# Why did EGM2008's quantization error come out 20x its design rms? Feed the
# same synthesis random coefficient errors of known size at a few latitudes,
# by degree band, and see where the excess comes from.
import sys, numpy as np
exec(open(sys.argv[2]).read().split("j = json.load")[0])   # constants, C/S, helpers
import json
rng = np.random.default_rng(1)
lat = np.array([0.0, 20.0, 40.0, 55.0]); lon = np.linspace(-180, 180, 721)[:-1]
hb = len(lat); phi = np.radians(lat); sp = np.sin(phi)
Nv = aw / np.sqrt(1 - e2 * sp ** 2); x = Nv * np.cos(phi); zc = Nv * (1 - e2) * sp
r = np.hypot(x, zc); tt = zc / r; uu = x / r
gamma = 9.7803253359 * (1 + 0.00193185265241 * sp ** 2) / np.sqrt(1 - e2 * sp ** 2)
q = a / r; lnq = np.log(q); lnu = np.log(uu); SCALE = 1e280
src = open(sys.argv[2]).read(); body = src[src.index('def synth'):src.index('t0 = time.time()')]
exec(body)
ncoef = (N + 1) * (N + 2) - (N + 1); d = 1e-3 * np.sqrt(12) / (a * np.sqrt(ncoef))
for lo, hi in ((2, 360), (361, 1000), (1001, 2190)):
    dC = np.zeros_like(C); dS = np.zeros_like(Sn)
    for n in range(lo, hi + 1):
        dC[n, :n + 1] = rng.uniform(-d / 2, d / 2, n + 1); dS[n, 1:n + 1] = rng.uniform(-d / 2, d / 2, n)
    dC[2191:] = 0
    err = synth(dC, dS)
    k = sum(2 * n + 1 for n in range(lo, hi + 1))
    expect = d / np.sqrt(12) * a * np.sqrt(k)
    print(f'degrees {lo}-{hi}: expected rms {expect*1000:.3f} mm; got by latitude', ' '.join(f'{la:g}deg:{np.sqrt(np.mean(e**2))*1000:.3f}' for la, e in zip(lat, err)), flush=True)
