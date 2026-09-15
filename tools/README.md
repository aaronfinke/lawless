# Diagnostics for neutron Laue scaling

Post-hoc analysis of a lawless unmerged output file (`UNMERGEDOUT`), used to
decide what the scale model is still missing.  All four need `gemmi` and
`numpy`, and an unmerged MTZ with a `LAMBDA` column.

The common idea: for every unique reflection with three or more observations,
compare each observation with the **leave-one-out** mean of its symmetry mates.
Resolution is fixed within a reflection, so a purely resolution-dependent error
cancels exactly and what survives is a genuine systematic.  The statistic is a
ratio of sums, `sum(I_obs) / sum(I_loo) - 1`, which is unbiased even when the
leave-one-out mean is noisy.

| script | question it answers |
|---|---|
| `twotheta_probe.py` | Where does the residual bias sit — in wavelength, resolution, or scattering angle? Prints a 2-theta x lambda table and both marginals. |
| `separability.py` | Which variable *owns* it? Fits a (lambda x d) grid of biases with a polynomial in lambda, in lambda and d, or in 2-theta, and reports the variance each explains. |
| `null_test.py` | Is the probe itself sound? Replaces observations by synthetic draws `I ~ N(mu_k, sigma_i^2)` keeping the real group structure, sigmas and labels, so no systematic exists by construction. Anything the probe reports here is an artefact. **Run this before believing any of the others.** |
| `twotheta_sim2.py` | What would an angle term buy? Fits `g(2theta)` against leave-one-out means and applies it once, then recomputes R-merge, R-pim, CC-half and chi-squared. |

Usage:

```
python twotheta_probe.py  unmerged.mtz "label"
python separability.py    "unmerged.mtz|label" ["unmerged2.mtz|label2" ...]
python null_test.py       unmerged.mtz [nreplicates]
python twotheta_sim2.py   unmerged.mtz "label" [chebyshev_degree]
```

Measured on CuZnSOD, dMPro-KB5 and hCAII-Cu (MaNDi, SNS), these found a
reproducible 2-theta systematic of about 9 % peak to peak which the
crystal-frame `SECONDARY` harmonics cannot absorb, and which a post-hoc
correction cannot fix — it has to be refined inside the scale model.  That is
the case for `SCALES ... TWOTHETA`.
