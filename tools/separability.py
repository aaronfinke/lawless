"""Which variable explains the residual bias: lambda, d, or 2-theta?

Build a (lambda x d) grid of leave-one-out biases, then ask how much of the
weighted variance across the grid is explained by
  (a) a smooth function of lambda alone,
  (b) an additive f(lambda) + g(d),
  (c) a smooth function of 2-theta alone (one parameter fewer than (b)).
"""
import sys
from collections import defaultdict
import numpy as np
import gemmi

def load(path):
    mtz = gemmi.read_mtz_file(path)
    cols = {c.label: i for i, c in enumerate(mtz.columns)}
    dat = np.array(mtz, copy=False).astype(np.float64)
    H, K, L = (dat[:, cols[c]].astype(int) for c in "HKL")
    I, S, LAM = (dat[:, cols[c]] for c in ("I", "SIGI", "LAMBDA"))
    asu, gops = gemmi.ReciprocalAsu(mtz.spacegroup), mtz.spacegroup.operations()
    groups = defaultdict(list)
    for i in range(len(H)):
        if S[i] > 0 and np.isfinite(I[i]) and LAM[i] > 0:
            groups[tuple(asu.to_asu((H[i], K[i], L[i]), gops)[0])].append(i)
    rec = []
    for idx in groups.values():
        if len(idx) < 3:
            continue
        idx = np.array(idx); Ii = I[idx]; n = len(idx)
        loo = (Ii.sum() - Ii) / (n - 1)
        sdloo = np.sqrt((S[idx]**2).sum() - S[idx]**2) / (n - 1)
        for j, k in enumerate(idx):
            dd = mtz.cell.calculate_d((H[k], K[k], L[k]))
            st = LAM[k] / (2*dd)
            if st < 1.0:
                rec.append((loo[j], sdloo[j], Ii[j], LAM[k], dd, np.degrees(2*np.arcsin(st))))
    r = np.array(rec)
    keep = r[:, 0] / r[:, 1] > 3.0
    return r[keep]

def analyse(path, label, nl=6, nd=6):
    r = load(path)
    loo, _, Iobs, lam, dsp, tt = r.T
    lq = np.quantile(lam, np.linspace(0, 1, nl+1))
    dq = np.quantile(dsp, np.linspace(0, 1, nd+1))
    bias, wt, Lc, Dc, Tc = [], [], [], [], []
    for a in range(nl):
        ml = (lam >= lq[a]) & ((lam < lq[a+1]) if a < nl-1 else (lam <= lq[nl]))
        for b in range(nd):
            m = ml & (dsp >= dq[b]) & ((dsp < dq[b+1]) if b < nd-1 else (dsp <= dq[nd]))
            if m.sum() < 60:
                continue
            bias.append(Iobs[m].sum()/loo[m].sum() - 1)
            wt.append(m.sum()); Lc.append(lam[m].mean()); Dc.append(dsp[m].mean()); Tc.append(tt[m].mean())
    bias, wt = np.array(bias), np.array(wt, float)
    Lc, Dc, Tc = map(np.array, (Lc, Dc, Tc))
    w = wt / wt.sum()
    def wrss(y, Xs):
        X = np.column_stack(Xs)
        W = np.sqrt(w)
        c = np.linalg.lstsq(X*W[:, None], y*W, rcond=None)[0]
        res = y - X@c
        return (w*res**2).sum(), X.shape[1]
    one = np.ones_like(bias)
    tot = (w*(bias - (w*bias).sum())**2).sum()
    def poly(x, k, xr):
        z = 2*(x - xr[0])/(xr[1]-xr[0]) - 1
        return [np.polynomial.chebyshev.chebval(z, np.eye(k+1)[i]) for i in range(1, k+1)]
    lamr, dr, ttr = (Lc.min(), Lc.max()), (Dc.min(), Dc.max()), (Tc.min(), Tc.max())
    models = {
        "lambda only (cubic)":        [one] + poly(Lc, 3, lamr),
        "lambda + d (both cubic)":    [one] + poly(Lc, 3, lamr) + poly(Dc, 3, dr),
        "2theta only (cubic)":        [one] + poly(Tc, 3, ttr),
        "2theta only (quintic)":      [one] + poly(Tc, 5, ttr),
        "2theta + lambda (cubic)":    [one] + poly(Tc, 3, ttr) + poly(Lc, 3, lamr),
    }
    print(f"\n=== {label} ===  {len(bias)} (lambda x d) cells, rms bias "
          f"{100*np.sqrt(tot):.2f}%")
    for name, Xs in models.items():
        rss, npar = wrss(bias, Xs)
        print(f"  {name:<26s} npar={npar:2d}  rms resid {100*np.sqrt(rss):5.2f}%   "
              f"var explained {100*(1-rss/tot):5.1f}%")

for spec in sys.argv[1:]:
    p, l = spec.split("|")
    analyse(p, l)
