"""Non-circular estimate of what a scattering-angle scale term would buy.

The correction is fitted against LEAVE-ONE-OUT means, so an observation never
influences the target it is fitted to, and applied in a single pass (no
alternation, hence no runaway feedback).  R-merge / R-pim are then recomputed.
"""
import sys
from collections import defaultdict
import numpy as np
import gemmi

path, label = sys.argv[1], sys.argv[2]
ncheb = int(sys.argv[3]) if len(sys.argv) > 3 else 4

mtz = gemmi.read_mtz_file(path)
cols = {c.label: i for i, c in enumerate(mtz.columns)}
dat = np.array(mtz, copy=False).astype(np.float64)
H, K, L = (dat[:, cols[c]].astype(int) for c in "HKL")
I, S, LAM = (dat[:, cols[c]] for c in ("I", "SIGI", "LAMBDA"))
asu, gops = gemmi.ReciprocalAsu(mtz.spacegroup), mtz.spacegroup.operations()

ks, Is, Ss, tts = [], [], [], []
lut = {}
for i in range(len(H)):
    if not (S[i] > 0 and np.isfinite(I[i]) and LAM[i] > 0):
        continue
    dd = mtz.cell.calculate_d((H[i], K[i], L[i]))
    st = LAM[i] / (2 * dd)
    if st >= 1.0:
        continue
    a = tuple(asu.to_asu((H[i], K[i], L[i]), gops)[0])
    ks.append(lut.setdefault(a, len(lut)))
    Is.append(I[i]); Ss.append(S[i]); tts.append(2 * np.arcsin(st))
k = np.array(ks); Io = np.array(Is); So = np.array(Ss); tt = np.array(tts)
ng = k.max() + 1
mult = np.bincount(k, minlength=ng)

def rfactors(Ic, Sc):
    w = 1.0 / Sc ** 2
    sw = np.bincount(k, weights=w, minlength=ng)
    m = np.bincount(k, weights=Ic * w, minlength=ng) / sw
    keep = mult[k] > 1
    den = m[k][keep].sum()
    n = mult[k][keep]
    rm = np.abs(Ic - m[k])[keep].sum() / den
    rp = (np.sqrt(1.0 / (n - 1)) * np.abs(Ic - m[k])[keep]).sum() / den
    return rm, rp

def cchalf(Ic, Sc, seed=0):
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(Ic))
    half = np.zeros(len(Ic), bool)
    seen = defaultdict(int)
    for i in perm:
        seen[k[i]] += 1
        half[i] = seen[k[i]] % 2 == 0
    out = []
    for sel in (half, ~half):
        w = 1.0/Sc[sel]**2
        sw = np.bincount(k[sel], weights=w, minlength=ng)
        m = np.zeros(ng); nz = sw > 0
        m[nz] = np.bincount(k[sel], weights=Ic[sel]*w, minlength=ng)[nz]/sw[nz]
        out.append((m, nz))
    (m1, n1), (m2, n2) = out
    both = n1 & n2
    return np.corrcoef(m1[both], m2[both])[0, 1], both.sum()

def chi2(Ic, Sc):
    w = 1.0/Sc**2
    sw = np.bincount(k, weights=w, minlength=ng)
    swi = np.bincount(k, weights=Ic*w, minlength=ng)
    l = (swi[k] - w*Ic)/(sw[k] - w)
    sl = np.sqrt(1.0/(sw[k] - w))
    ok2 = mult[k] > 1
    return np.mean(((Ic - l)**2/(Sc**2 + sl**2))[ok2])

def loo(Ic, Sc):
    w = 1.0 / Sc ** 2
    sw = np.bincount(k, weights=w, minlength=ng)
    swi = np.bincount(k, weights=Ic * w, minlength=ng)
    return (swi[k] - w * Ic) / (sw[k] - w)

c2 = np.cos(tt)
x = 2 * (c2 - c2.min()) / (c2.max() - c2.min()) - 1
B = np.column_stack([np.polynomial.chebyshev.chebval(x, np.eye(ncheb + 1)[j])
                     for j in range(1, ncheb + 1)])

r0 = rfactors(Io, So)
ref = loo(Io, So)
good = (mult[k] > 2) & (ref > 3 * So)            # selection on the LOO value only
y = np.log(Io[good] / ref[good])
ww = np.sqrt(1.0 / (So[good] / ref[good]) ** 2)  # weight by loo I/sigma
a = np.linalg.lstsq(B[good] * ww[:, None], y * ww, rcond=None)[0]
lg = B @ a
lg -= np.average(lg[good], weights=ww ** 2)
g = np.exp(lg)
r1 = rfactors(Io / g, So / g)

print(f"\n=== {label} ===  degree {ncheb} in cos2theta, {ncheb} free parameters")
print(f"  R-merge  {r0[0]:.4f} -> {r1[0]:.4f}  ({100*(r1[0]-r0[0])/r0[0]:+.2f} %)")
print(f"  R-pim    {r0[1]:.4f} -> {r1[1]:.4f}  ({100*(r1[1]-r0[1])/r0[1]:+.2f} %)")
c0, nn = cchalf(Io, So); c1, _ = cchalf(Io/g, So/g)
print(f"  CC-half  {c0:.4f} -> {c1:.4f}   ({nn} unique)")
print(f"  chi^2    {chi2(Io,So):.3f} -> {chi2(Io/g,So/g):.3f}")
deg = np.degrees(tt)
q = np.quantile(deg, np.linspace(0, 1, 7))
print("  fitted g(2theta):  " + "  ".join(
    f"{q[r]:.0f}-{q[r+1]:.0f}deg {g[(deg>=q[r])&(deg<=q[r+1])].mean():.3f}" for r in range(6)))
