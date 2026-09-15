"""Is the residual bias a function of scattering angle 2-theta?

2-theta = 2 asin(lambda / 2d) is the one variable that couples wavelength and
resolution.  If the (d x lambda) interaction seen in bias_probe2 collapses onto
2-theta alone, the missing term in the scale model is angle-dependent
(absorption path / detector obliquity / panel), not wavelength-dependent.
"""
import sys
from collections import defaultdict
import numpy as np
import gemmi

path, label = sys.argv[1], sys.argv[2]
mtz = gemmi.read_mtz_file(path)
cols = {c.label: i for i, c in enumerate(mtz.columns)}
d = np.array(mtz, copy=False).astype(np.float64)
H, K, L = (d[:, cols[c]].astype(int) for c in "HKL")
I, S, LAM = (d[:, cols[c]] for c in ("I", "SIGI", "LAMBDA"))
sg, cell = mtz.spacegroup, mtz.cell
asu, gops = gemmi.ReciprocalAsu(sg), sg.operations()

groups = defaultdict(list)
for i in range(len(H)):
    if S[i] > 0 and np.isfinite(I[i]) and LAM[i] > 0:
        groups[tuple(asu.to_asu((H[i], K[i], L[i]), gops)[0])].append(i)

rec = []
for idx in groups.values():
    if len(idx) < 3:
        continue
    idx = np.array(idx)
    Ii = I[idx]
    n = len(idx)
    loo = (Ii.sum() - Ii) / (n - 1)
    sdloo = np.sqrt((S[idx] ** 2).sum() - S[idx] ** 2) / (n - 1)
    for j, k in enumerate(idx):
        dd = cell.calculate_d((H[k], K[k], L[k]))
        st = LAM[k] / (2 * dd)
        if st >= 1.0:
            continue
        rec.append((loo[j], sdloo[j], Ii[j], LAM[k], dd, np.degrees(2 * np.arcsin(st))))

rec = np.array(rec)
loo, sdloo, Iobs, lam, dsp, tt = rec.T
keep = loo / sdloo > 3.0
loo, Iobs, lam, dsp, tt = loo[keep], Iobs[keep], lam[keep], dsp[keep], tt[keep]

print(f"\n=== {label} ===   n={len(loo)}   2theta {tt.min():.1f}-{tt.max():.1f} deg")

def table(var, name, nr, fmt):
    lq = np.quantile(lam, np.linspace(0, 1, 6))
    vq = np.quantile(var, np.linspace(0, 1, nr + 1))
    print(f"\n  bias %  [rows {name}, cols lambda]  " +
          "".join(f"{lq[c]:.2f}-{lq[c+1]:.2f} " for c in range(5)))
    for r in range(nr):
        m = (var >= vq[r]) & ((var < vq[r+1]) if r < nr - 1 else (var <= vq[nr]))
        line = f"  {vq[r]:{fmt}}-{vq[r+1]:<{fmt[:-1]}}"
        for c in range(5):
            mc = m & (lam >= lq[c]) & ((lam < lq[c+1]) if c < 4 else (lam <= lq[5]))
            line += f"{100*(Iobs[mc].sum()/loo[mc].sum()-1):>8.2f} " if mc.sum() > 100 else "    -    "
        print(line + f"  n={m.sum()}")

table(tt, "2theta (deg)", 5, "6.1f")

print("\n  marginal bias vs 2theta alone")
q = np.quantile(tt, np.linspace(0, 1, 9))
for r in range(8):
    m = (tt >= q[r]) & ((tt < q[r+1]) if r < 7 else (tt <= q[8]))
    print(f"    {q[r]:6.1f}-{q[r+1]:6.1f} deg : {100*(Iobs[m].sum()/loo[m].sum()-1):+7.2f} %   n={m.sum()}")

print("\n  marginal bias vs lambda alone")
q = np.quantile(lam, np.linspace(0, 1, 9))
for r in range(8):
    m = (lam >= q[r]) & ((lam < q[r+1]) if r < 7 else (lam <= q[8]))
    print(f"    {q[r]:6.2f}-{q[r+1]:6.2f} A   : {100*(Iobs[m].sum()/loo[m].sum()-1):+7.2f} %   n={m.sum()}")
