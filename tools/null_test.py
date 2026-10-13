"""Null test for the leave-one-out bias probe.

Replace every observation by a synthetic draw I_i ~ N(m_k, sigma_i^2), keeping the
real group structure, the real sigmas, and the real (lambda, d, 2theta) labels.
There is now NO systematic by construction.  If the probe still reports a 2theta
trend, the trend is an artefact of the probe, not of the data.
"""
import sys
import numpy as np
import gemmi

path = sys.argv[1]
nrep = int(sys.argv[2]) if len(sys.argv) > 2 else 3

mtz = gemmi.read_mtz_file(path)
cols = {c.label: i for i, c in enumerate(mtz.columns)}
dat = np.array(mtz, copy=False).astype(np.float64)
H, K, L = (dat[:, cols[c]].astype(int) for c in "HKL")
I, S, LAM = (dat[:, cols[c]] for c in ("I", "SIGI", "LAMBDA"))
asu, gops = gemmi.ReciprocalAsu(mtz.spacegroup), mtz.spacegroup.operations()

ks, Is, Ss, tts, lut = [], [], [], [], {}
for i in range(len(H)):
    if not (S[i] > 0 and np.isfinite(I[i]) and LAM[i] > 0):
        continue
    dd = mtz.cell.calculate_d((H[i], K[i], L[i]))
    st = LAM[i] / (2 * dd)
    if st >= 1.0:
        continue
    ks.append(lut.setdefault(tuple(asu.to_asu((H[i], K[i], L[i]), gops)[0]), len(lut)))
    Is.append(I[i]); Ss.append(S[i]); tts.append(np.degrees(2 * np.arcsin(st)))
k = np.array(ks); Io = np.array(Is); So = np.array(Ss); tt = np.array(tts)
ng = k.max() + 1
mult = np.bincount(k, minlength=ng)

def probe(Ivals):
    """the same leave-one-out ratio-of-sums bias, marginal in 2theta"""
    n = mult[k].astype(float)
    tot = np.bincount(k, weights=Ivals, minlength=ng)[k]
    loo = (tot - Ivals) / (n - 1)
    s2 = np.bincount(k, weights=So ** 2, minlength=ng)[k]
    sdloo = np.sqrt(s2 - So ** 2) / (n - 1)
    m = (n > 2) & (sdloo > 0) & (loo > 3 * sdloo)
    q = np.quantile(tt[m], np.linspace(0, 1, 9))
    out = []
    for r in range(8):
        b = m & (tt >= q[r]) & (tt <= q[r + 1])
        out.append(100 * (Ivals[b].sum() / loo[b].sum() - 1))
    return np.array(out), q

real, q = probe(Io)
mean = np.bincount(k, weights=Io, minlength=ng) / mult
rng = np.random.default_rng(0)
sims = np.array([probe(mean[k] + rng.normal(0, So))[0] for _ in range(nrep)])

print(f"\n{path}")
print(f"{'2theta bin':>20s}  {'real':>8s}  {'null mean':>10s}  {'null sd':>8s}  {'excess':>8s}")
for r in range(8):
    print(f"  {q[r]:7.1f}-{q[r+1]:7.1f}  {real[r]:8.2f}  {sims[:,r].mean():10.2f}  "
          f"{sims[:,r].std():8.2f}  {real[r]-sims[:,r].mean():8.2f}")
