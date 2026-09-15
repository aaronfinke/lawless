# aimless — CCP4 Data Scaling Program

## Overview

`aimless` scales and merges X-ray (and now Laue) crystallographic reflection data from unmerged MTZ files. Written in C++ using CCP4/Clipper libraries.

**Build:**
```bash
cmake --preset default        # configure (uses CMakePresets.json)
make -C build -j4             # compile
```

CMakePresets.json in the project root defines all library paths. Headers come from `~/ccp4-8.0-src/checkout`; dylibs from `/Applications/ccp4-9/lib`; `libcctbx.dylib` from `~/miniforge3/envs/cctbx/lib`. A POST_BUILD rule patches Jenkins-originated rpath entries automatically.

In VSCode: `Cmd+Shift+B` runs the default build task.

---

## Scale model

Each reflection's scale factor is:

```
g = ps × bs × ss × ds × ws
```

| Factor | Symbol | Description |
|--------|--------|-------------|
| Primary scale | `ps` | Per-run/batch multiplicative scale |
| B-factor | `bs` | Isotropic B via `exp(-B·s²)` |
| Secondary beam | `ss` | Angular correction |
| Detector/tiling | `ds` | Detector gain correction |
| **Wavelength** | `ws` | Chebyshev *or* GP normalisation, optionally times a per-run residual (Laue only) |

`ws = 1.0` when Laue mode is inactive (monochromatic data; backward-compatible).

---

## Chebyshev wavelength normalisation

### Parameterisation — log space

Coefficients represent **h(λ) = log f(λ)** fitted as a Chebyshev polynomial:

```
z(λ) = (2λ − λmin − λmax) / (λmax − λmin)   [maps [λmin,λmax] → [−1,1]]
h(λ) = Σ_{k=0}^{N} a_k T_k(z)               [Clenshaw's algorithm on LOG scale]
ws(λ) = exp(h(λ) − h(λ_ref))                 [always positive by construction]
```

Using the log of the scale polynomial (rather than the scale directly) guarantees `ws > 0` for any coefficient values, preventing the BFGS from stepping into singular or negative-scale regions.

Derivative for least-squares refinement:

```
dws/da_k = ws · (T_k(z) − T_k(z_ref))
```

Initial coefficients `[0, 0, 0, ...]` give `ws = 1` everywhere (flat spectrum).

Up to 5 wavelength ranges supported, each with its own polynomial degree and coefficient set. Ranges are non-overlapping; reflections outside all ranges get `ws = 1.0`.

### Class: `WavelengthChebyshevScale` (`scaletypes.hh/.cpp`)

```cpp
WavelengthChebyshevScale(
    const std::vector<WavelengthRange>& Ranges,
    const double& LambdaRef);

struct WavelengthRange {
    double lam_min, lam_max;
    int degree;
    int offset;  // index into flat coefficients vector (set by constructor)
};

double Scale(const double& lambda) const;          // returns ws = exp(h(λ) - h_ref)
double ScaleDeriv(const double& lambda,
                  std::vector<double>& dgdp) const; // fills gradient, returns ws
std::string PrintNormalization(int npoints=12) const; // log table for output file
```

`f_ref` in the private members stores `h(λ_ref)` (the log-polynomial value at the reference wavelength), not f(λ_ref) as the name might suggest.

### Integration in `ScaleModel` (`scalemodel.hh/.cpp`)

- Single global `WavelengthChebyshevScale wavelength_scale` per `ScaleModel`
- Parameters appended after detector block in the flat parameter vector (`idxwavscale`)
- `GetParameterType()` returns `ScaleModel::WAVELENGTH` for these indices
- `HasWavelengthScale()` — predicate used to trigger the Laue pre-pass
- `SetWavelengthOnlyMode(bool)` — when true, `ScaleFactorDeriv` zeros all non-WAVELENGTH derivatives so BFGS only moves Chebyshev coefficients
- `GetLargeShift()` returns 0.5 for WAVELENGTH parameters (limits BFGS step size)
- `IsRefinable()` includes `nwavscale > 0` so Laue-only mode is not treated as non-refinable

---

## Gaussian-process (GPR) wavelength normalisation

Non-parametric alternative to the Chebyshev model. Selected with `LAUE NORMGPR`.
**Mutually exclusive with `NORMCHEBYSHEV`** — the parser raises a fatal error if
both are given (in either order). Use one method or the other.

### Key difference from Chebyshev: it is NOT refined in BFGS

The GP is fitted **once** in a dedicated pre-pass directly from the data, then
applied as a **fixed multiplicative correction** (a precomputed lookup table).
It contributes **zero parameters** to the scale-model parameter vector and has
**zero derivatives** in `ScaleFactorDeriv`. So none of the parameter-management
machinery (`CountParameters`, `GetParameters/SetParameters`, `GetParameterType`,
`GetLargeShift`) is touched — much simpler than the Chebyshev integration.

### Algorithm (`WavelengthGPRScale::Fit`, `scaletypes.cpp`)

1. **Empirical response** (built in `aimless.cpp`, not via BFGS): for each
   reflection with ≥2 accepted observations passing the I/σ cut, form the
   **leave-one-out** ratio `ρ = I_obs / <I>_other-mates` at each observation's
   wavelength, with weight `u = <I>_other-mates / σ²`. Leave-one-out avoids the
   self-bias that a plain mean would introduce at low multiplicity. Weak and
   negative intensities are kept — the ratio estimator handles them, and
   cutting them biases the response at the spectrum edges.
   The whole sample-building step is **iterated 3×**: each pass divides the
   mates by the current `ws(λ)` before averaging, so the mates' mean is free of
   wavelength response. Without this the raw ratio measures the response only
   relative to the (wavelength-dependent) mixture of wavelengths at which the
   mates happen to have been measured. The iteration converges on the same
   consistency condition the Chebyshev refinement reaches by least squares
   (typically converged by pass 2).
2. **Bin** the sorted `(λ, ρ, u)` samples into **equal-count** bins, closing a
   bin at whichever comes first: the nominal population `N/nbins`, or a cap of
   `4·(log-λ span)/nbins` on the bin's width in `log λ`. Equal count alone
   gives every training point comparable precision (uniform-width bins do not —
   a Laue dataset has orders of magnitude more observations at the peak of the
   spectrum than in its tails); the width cap stops the sparse tail collapsing
   into one very wide bin whose centroid sits far from the end of the range.
   Each bin gives one training point by the linear-space **ratio estimator**
   `r = Σuρ/Σu` with the sandwich variance `Σu²(ρ−r)²/(Σu)²`; the target is
   `h = log r`, `var(h) = var(r)/r²`. A mean of per-observation *log*-ratios
   (the original estimator) is badly biased and heavy-tailed for weak data.
3. **Fit** a GP to `h_b` in log-intensity space with a squared-exponential
   (default) or Matérn-3/2 kernel, **a constant mean**, and the kernel
   evaluated in the **warped coordinate `x = log λ`**. Length scale `ℓ`,
   signal sdev `σ_f` and a **noise-inflation factor `α`** (multiplying the bin
   variances) are chosen *together* by maximising the **log marginal
   likelihood** over a grid. Solved via Eigen `LLT` (Cholesky) of
   `σ_f²K + diag(α·var) + jitter`.
4. **Lookup table**: evaluate the posterior mean `g(λ)` on a dense uniform grid,
   **held constant outside the outermost training points** so the fit can never
   run away where there is no data. The posterior SD is still evaluated at the
   true λ, so the widening band keeps showing that the extremes are unsupported.

`Scale(λ) = exp(g(λ) − g(λ_ref))` — log space guarantees `ws > 0`; reference
normalisation makes `ws(λ_ref) = 1`. Returns `1.0` outside `[λmin, λmax]` or if
the fit was skipped (too few samples/bins).

#### Why these choices (measured, not assumed)

The first implementation produced strongly oscillatory `w(λ)`. Diagnosis on
`test_data/ca_thio`: the marginal likelihood pinned `ℓ` at the **lower** end of
its search grid on nearly every trial, i.e. the GP was interpolating bin means.
Two causes, both fixed above:

- **Noise underestimated.** `SEM² = var/count` ignored that `1/σ²` weights make
  the *effective* sample size several times smaller than the count (Kish
  `n_eff = (Σw)²/Σw²` was as low as 5 for a 151-observation bin) — a noise
  variance too small by up to ~25×.
- **A single stationary kernel in λ cannot fit both ends.** The spectrum rises
  steeply over a narrow interval at short λ and is broad and flat at long λ;
  short `ℓ` was the only way to follow the rise, and it made everything else
  wiggle. `x = log λ` makes the response near-stationary.

Held-out cross-validation (fit on a random half of the observations, score
against bins from the other half) gives **χ²/n = 4.34 for the original fit vs
2.22 for the current one** — the oscillations did not reproduce on independent
data. Note that the original fit gives a slightly *lower* Rmerge (0.218 vs
0.226), which is what overfitting looks like: R-factors reward a curve that
absorbs random variation. On the same test the current GP matches the Chebyshev
model (Rmerge 0.226 vs 0.224, Rmeas 0.276 both, `<I/σ>` 14.6 vs 14.5) while
staying smooth and flat where the Chebyshev diverges at long λ.

### Integration in `ScaleModel`

- Members `WavelengthGPRScale gpr_scale`, `GPRControl gpr_control`,
  `bool gpr_requested`, `double gpr_lambda_ref`, and for the per-run option
  `std::vector<WavelengthGPRScale> gpr_run_scales`, `bool gpr_perrun`.
- `HasGPRWavelengthScale()` — triggers the GP pre-pass (parallel to
  `HasWavelengthScale()` for Chebyshev).
- `FitGPRWavelength(lambdas, logratios, weights, output)` — driver called from
  the pre-pass; forwards to `gpr_scale.Fit`.
- In `ScaleFactor`/`ScaleFactorDeriv`: `ws *= gpr_scale.Scale(obs.lambda())`,
  then `*= gpr_run_scales[jscale].Scale(obs.lambda())` when the per-run residual
  is active (`jscale` is the run index — see `GPRRunActive`).
  In the deriv path the fixed factor `wsgpr` also scales the Chebyshev
  derivative block (defensive: keeps the two correct if ever combined, though
  the parser forbids that — they are mutually exclusive keywords).
- `haveWavelength` in `init` is true when **either** Chebyshev ranges or GPR is
  requested, so the constant-primary fallback (LAMBDAONLY without
  `SCALES CONSTANT`) works for GPR too.

### Per-run residual normalisation (`LAUE NORMGPRPERRUN`)

A single global `w(λ)` is a compromise if the incident spectrum drifts between
runs. `NORMGPRPERRUN` fits each run an **additional residual** curve, so the
total correction is

```
ws(λ, run) = w_global(λ) × w_run(λ)
```

Fitting a residual rather than a full per-run spectrum is what keeps this well
conditioned: the global curve carries the (large) spectral shape using all the
data, and each run only has to explain its departure from it. A run with no
departure fits flat at 1.0 and changes nothing; a run with too few samples
degrades to the global curve rather than to noise.

The residual samples are the same leave-one-out ratios as the global fit, with
one change: **both the observation and its mates are first divided by the total
correction already in force** (`ScaleModel::GPRWavelengthScaleAt(irun, λ)`).
Samples are assigned to the run of the *observation*, while the mates may come
from any run — that is what makes the residual relative to the common
consensus rather than to the run itself. Two iterations (the global fit uses
three).

Because the GP contributes **no parameters** to the refinement vector, this
needed no parameter-management changes at all — only the vector of
`WavelengthGPRScale` indexed by run and the extra factor in the two scale
paths.

- `HasPerRunGPRWavelength()`, `GPRRunActive(irun)`, `NumberGPRRuns()`
- `FitGPRWavelengthRun(irun, lambdas, ratios, weights, output)` — per-run driver
- `PrintGPRPerRunNormalization(output)` — logs each run's residual table
- The per-run loop lives in `aimless.cpp`, after the global GP loop, and only
  runs if the global fit succeeded.

**Diagnosing a null result:** compare the per-run `sigma_f` in the log against
the global one. A per-run `sigma_f` two or more orders of magnitude smaller
means the GP is reporting that there is no per-run drift to correct — which is
a real answer, not a failed fit. Note also that this corrects *means*: a
per-run wavelength effect that is an increase in **scatter** rather than a
systematic shift will fit flat, and no wavelength normalisation of any kind can
absorb it (only the SD correction can).

### Log & XML output

`WavelengthGPRScale::PrintNormalization` (via
`ScaleModel::PrintGPRWavelengthNormalization`) writes to the log: the reference
wavelength, range, kernel, length scale (in `log λ` units **and** its Å
equivalent at `λ_ref`), `σ_f`, noise-inflation factor, training-bin count, a
`w(λ)` sample table, and an ASCII line plot of `w(λ)`
(`WavelengthGPRScale::AsciiPlot`, reference wavelength drawn as a `:` column).
Each of the 3 fit passes logs its own one-line summary, plus a warning if the
length scale lands on either end of its search grid (a sign the fit is
following noise, or is featureless).

`WavelengthGPRScale::asXML()` (via `ScaleModel::GPRWavelengthNormalizationXML()`,
emitted from `aimless.cpp` after the log table) writes a
`<WavelengthNormalisationGPR>` block to XMLOUT: `<ReferenceWavelength>`,
`<LambdaMin>`/`<LambdaMax>`, `<Kernel>`, `<LengthScale>`, `<SigmaF>`,
`<NoiseInflation>`, `<TrainingBins>`, and a `<Normalisation>` table of 21
`<point>` (`<lambda>`,`<w>`,`<uncertainty>`) samples. Mirrors the Chebyshev
`<WavelengthNormalisation>` block on the `lawless` branch.

### LAMBDANORM gnuplot file + uncertainties

The GP fit also produces a **posterior SD** (`grid_sd`, computed in `Fit()` from
the Cholesky factor: `var = k(λ,λ) − v·v`, `v = L⁻¹k_*`).
`WavelengthGPRScale::Uncertainty(λ)` interpolates it; in log space this is the
relative (fractional) uncertainty of `ws(λ)`. It appears as `<uncertainty>` in
the XML and as column 3 of the LAMBDANORM file.

`WavelengthGPRScale::GnuplotScript(title, version)` (via
`ScaleModel::GPRWavelengthGnuplot()`) returns a **self-contained gnuplot script**
written by `aimless.cpp` to a file named **`LAMBDANORM`** (in the GP pre-pass).
It has a header comment (program/version from `version.hh`, run title, and the
`gnuplot -p LAMBDANORM` open instruction), an inline `$LAMBDANORM` datablock
(columns: `λ  w  rel_uncertainty  w_lo  w_hi`), a `$LAMBDABINS` datablock with
the **binned observations the GP was fitted to** (`λ  w_bin  σ(w_bin)`), and a
`plot` of the `w(λ)` line over a 1σ `filledcurves` band with the bins as
error-bar points. The bins are the diagnostic that matters: the curve should
follow them without chasing individual points.
Gnuplot (≥5.0) is preferred over the deprecated
loggraph plot files (ROGUES/SCALES/ANOMPLOT/CORRELPLOT). GPR-only; the
Chebyshev model writes no LAMBDANORM.

Advantage over Chebyshev: adapts to variable data density via the
heteroscedastic noise model and avoids polynomial runaway at sparse wavelength
extremes (outside the range it returns a flat `ws = 1`, never a divergent tail).

Eigen3 lives at `${SRC}/eigen-3.4.0` and is added to `include_directories` in
`CMakeLists.txt`; the dependency is confined to `scaletypes.cpp`.

---

## Secondary-beam / absorption geometry for Laue

The secondary beam is `ŝ_out = ŝ_in + λ·q`, so the **wavelength sets the
scattering angle**. `Batch::HtoSr0` used to form the diffraction vector from
`batchinfo.alambd`, the batch wavelength — exact for monochromatic data, where
every observation in the batch shares it, and wrong for Laue. Over a 2.10–3.90 Å
band that misdirects the secondary beam by a **median 6.4°**, a quarter of
observations by more than 10°, worst case 18°. The azimuth about the beam is
unaffected; the whole error is in 2θ and is systematic in λ.

`HtoSr0`, `CalcSecondaryBeam` and `CalcSecondaryBeamPolar` now take an optional
wavelength, defaulting to `-1` meaning "use the batch wavelength", so every
existing call is unchanged. The single caller, `CalcSecondaryBeams`, already
holds the observation and passes `observation::lambda()` — which itself falls
back to the batch wavelength when the file has no LAMBDA column, so
monochromatic behaviour is preserved twice over.

Each of those three functions has exactly one call site; the values they produce
are consumed only by `obs.GetS2()` → the secondary/absorption scale in
`scalemodel.cpp`, and `obs.GetS()` → the axes of the ROGUES plot.

### What the batch header needs before any of this works

`Batch::init` rejects orientation data if `det(U) <= 0.9` or `U == I`, and
`ScaleModel` then refuses `SECONDARY`/`ABSORPTION` with "can't use Secondary
Scale unless we have valid orientation information". To supply it:

- `UMAT` is **column-major** (`MVutil::SetCMat33` reads `m[0],m[3],m[6]` as the
  first row), and must satisfy `U·B = UB` with `B` exactly as `Scell::Bmat`
  builds it from the **symmetry-constrained** cell — pointless constrains batch
  cells on output, so a `U` paired with an unconstrained cell will not survive.
- `SOURCE` is **anti-parallel to the beam**, `E1` is the principal rotation
  axis. `MakePermutationMatrix(SOURCE, E1)` builds the Cambridge frame from
  those two, putting z exactly on `E1` and x in the beam/`E1` plane, so the
  header frame need not be canonical — it only has to be self-consistent and
  to name the beam and the axis correctly.
- **A `ROT` column is mandatory.** `mtz_unmerge_io.cpp` has
  `if (col_sel.col_Rot < 0) {phi = batch;}` — with no ROT column aimless uses
  the *batch number* as φ and back-rotates the source vector by hundreds of
  degrees. For stills whose exposure-to-exposure rotation is carried in `U` and
  `SOURCE`, write `ROT = 0`.

---

## Wavelength term in the SD correction (`SDCORRECTION SDLAMBDA`)

`SDCORRECTION` fits `sigma' = SdFac*sqrt(sigma^2 + SdB*I + (SdAdd*I)^2)` — three
terms, all functions of **intensity**. Laue data can carry a variance component
that depends on **wavelength** instead, which none of them can represent.

It is a variance, not a bias, so no wavelength normalisation can absorb it
either: `LAUE NORMGPRPERRUN` correctly fits a flat residual against it, because
a scale model corrects means. Giving the intensity-only model its own copy per
run does not help either — `SDCORRECTION REFINE INDIVIDUAL` sends SdFac to
0.46–1.77 and SdAdd to 0.10–0.61 in opposition (they trade off directly, since
SdFac multiplies the SdAdd term) and trebles the pack-to-pack chi^2 spread. The
missing ingredient is the *shape*, not the granularity.

```
SDCORRECTION ... SDLAMBDA [<scale>]     ->   sigma' *= (lambda/lambda_ref)^SdLam
```

One exponent per run, from `SDmodel::FitLambdaSDcorrection`. chi^2 is a variance
ratio, so scaling sigma by F divides it by F^2; with `chi^2 ~ (lambda/lref)^m`
the exponent is `m/2`, and the constant part is left to SdFac. `lambda_ref` is
`LAUE NORMLAMREF`. The optional `<scale>` multiplies the fitted exponent and
defaults to 1 — it is a manual damping knob, not normally needed.

### Four things the fit has to get right

These were all found the hard way; three of them produce a *plausible-looking*
wrong answer rather than an obvious failure.

1. **Scaled intensities.** Use `obs.kI()` / `obs.ksigI()` throughout, including
   the leave-one-out mean. Mixing scaled and unscaled is wrong anywhere, but for
   Laue it is wrong in a specifically wavelength-correlated way, because the
   scale is itself a strong function of wavelength — it manufactures exactly
   the trend the term then "corrects".
2. **Fit within resolution bins.** lambda and resolution are correlated in Laue
   data (`lambda = 2d sin(theta)`) and chi^2 has a strong resolution trend of
   its own, so a raw slope against lambda silently absorbs part of it. The fit
   bins by (resolution x lambda) and pools residuals centred within each
   resolution bin.
3. **Fit late.** It runs after the merge-stage outlier rejection in
   `aimless.cpp`, not inside `AnalyseSD`. The residual before final scaling and
   rejection is not the residual in the delivered data.
4. **Leave-one-out.** The deviation is of an observation from the mean of its
   *mates*, not from a mean it helped define.

With 1 and 3 wrong the correction was not self-limiting: on a dataset with a
flat residual it fitted exponents of -0.3 to -0.65 and *created* a wavelength
trend.

### Checking it on new data

It should be a no-op when there is nothing to correct. Verified on two datasets:

| dataset | residual | fitted exponents | outcome |
|---|---|---|---|
| dMPro-KB5 (C2, 15 packs, mult 3.4) | flat | -0.04 to +0.11 | chi^2 vs lambda unchanged to 2 dp |
| CuZnSOD (P6522, 10 packs, mult 11.5) | +0.23 across the band | 0.12 to 0.35 | chi^2 vs lambda flattened to -0.02 |

Merging statistics barely move either way — sigma enters the merged mean only
through the weights. The value is error estimates that mean what they say, for
refinement weights and resolution cuts, not better R-factors.

**Before trusting it on a new dataset, check the residual is a variance and not
a bias**: plot mean `I/<mates>` against lambda alongside chi^2 against lambda.
If the *means* trend, the scaling is at fault and inflating sigma would hide it.
Only if the means are flat and chi^2 trends does the error model own it.

A single power law cannot represent a non-monotone wavelength dependence (a
detector edge, a bandpass edge). If one turns up, `WavelengthGPRScale` is the
natural generalisation — the same fit applied to log sigma instead of log I.

---

## Workflow for Laue data

### Wavelength pre-normalisation pass

Before the normal first-round and main scaling, aimless runs a dedicated wavelength-only pass:

1. Reflections selected by the same I/σ criterion as first-round scaling
2. `SetWavelengthOnlyMode(true)` — only Chebyshev coefficients are refined (BFGS sees zero gradient for ps, bs, ss, ds)
3. 5 cycles of BFGS to establish the wavelength normalization curve
4. `PrintWavelengthNormalization()` writes the curve to the log file
5. `SetWavelengthOnlyMode(false)` — full joint refinement resumes normally

**Why this ordering matters:** primary scale (per-batch) and wavelength scale are strongly correlated for Laue data (each batch sees the full wavelength spectrum). Refining wavelength first decorrelates them, giving better starting values for the subsequent joint refinement.

### Wavelength normalisation log output

The log shows, for each range:
- Log-coefficients `a_k`
- Table of `w(λ)` at 12 sample points across the range
- Reference wavelength marked with `<- ref`
- An ASCII line plot of `w(λ)` (via `WavelengthChebyshevScale::AsciiPlot`),
  with the reference wavelength drawn as a `:` column

### Wavelength normalisation XML output

`WavelengthChebyshevScale::asXML()` (surfaced through
`ScaleModel::WavelengthNormalizationXML()`) emits a `<WavelengthNormalisation>`
block to the XMLOUT file, written from `aimless.cpp` right after the log table.
It contains `<ReferenceWavelength>` and, per `<Range>`, the `<LambdaMin>`,
`<LambdaMax>`, `<Degree>`, `<LogCoefficients>`, and a `<Normalisation>` table of
21 `<point>` (`<lambda>`,`<w>`) samples.

### Wavelength normalisation gnuplot file (LAMBDANORM)

`WavelengthChebyshevScale::GnuplotScript(title, version)` (via
`ScaleModel::WavelengthGnuplot()`) returns a self-contained gnuplot script
written by `aimless.cpp` to a file named **`LAMBDANORM`** in the wavelength
pre-pass. Header comment (program/version, run title, `gnuplot -p LAMBDANORM`
hint), explicit `set xrange`/`set yrange`, then per-range datablocks and a
`plot`. Each range is drawn **in full across the whole domain in its own
colour**, using that range's own polynomial (`evalRange`): **solid** within the
range's interval (`$LNs<ir>`, the used part) and **dashed** where the polynomial
is extrapolated outside its interval (`$LNd<ir>`, not used). `yrange` is taken
from the used portions only, so divergent extrapolations run off-screen and are
clipped; non-finite extrapolated points are skipped. A single range collapses to
one solid `w(λ)` line.
**No uncertainty band** — the Chebyshev fit carries no posterior covariance
(unlike the GP, whose LAMBDANORM has a 1σ band). Both methods write to the same
`LAMBDANORM` filename, but `NORMCHEBYSHEV` and `NORMGPR` are mutually exclusive
so only one path runs.

### Complete scaling workflow (Laue)

```
Initial scaling
  ↓
Wavelength pre-normalisation (5 cycles, WAVELENGTH params only)
  ↓
First-round scaling (CONSTANT scales, secondary off, symmetric tiles)
  ↓
Outlier rejection
  ↓
Main scaling (full model, secondary on, asymmetric tiles)
  ↓
Merge and output
```

With `LAMBDAONLY` the workflow collapses to just the wavelength pass:

```
Wavelength normalisation (REFINE CYCLES, WAVELENGTH params only)
  ↓
Apply ws-only scales (no outlier rejection)
  ↓
Merge and output (unmerged MTZ)
```

---

## Using Laue wavelength scaling

### Keywords

```
LAUE NORMCHEBYSHEV <degree> <lam_min> <lam_max>
LAUE NORMLAMREF <lambda_ref>
```

- `NORMCHEBYSHEV` — add a wavelength range (repeat for each range, up to 5)
- `NORMLAMREF` — reference wavelength; `ws(λ_ref) = 1.0` by definition

**Note:** both sub-keywords start with "NORM", but `keyIs()` only compares 4 characters. The parser disambiguates on character 5 (`C` vs `L`).

#### GPR (Gaussian-process) alternative

```
LAUE NORMGPR <lam_min> <lam_max>   # enable GP normalisation over this range
LAUE NORMGPRLENGTH <lengthscale>   # optional: fix GP length scale (Å); default auto
LAUE NORMGPRBINS <nbins>           # optional: number of training bins; default 50
LAUE NORMGPRMATERN                 # optional: Matérn-3/2 kernel; default squared-exp
LAUE NORMGPRPERRUN                 # optional: residual w(λ) per run, on top of the global curve
LAUE NORMLAMREF <lambda_ref>       # shared with Chebyshev
```

`NORMGPRBINS` is the *nominal* bin count — the width cap can add a few bins in
sparse regions, so the reported "training bins" may slightly exceed it. Fewer
bins is the lever if a fit still looks too free.
`NORMGPRLENGTH` is given in Å but the kernel works in `log λ`, so the value is
converted as `ℓ_x = ℓ / λ_ref`; the log reports both forms.

The `NORMGPR*` sub-keywords share the 4-char key `NORM` with `keyIs()`, so the
parser disambiguates on character 5 (`G`) and then on the **full token string**
(`NORMGPR` vs `NORMGPRLENGTH` vs `NORMGPRBINS` vs `NORMGPRMATERN` vs
`NORMGPRPERRUN`), since
`keyIs("NORMGPR")` would also match the longer variants. One range only (the GP
spans the whole `[lam_min, lam_max]`); no `degree`. Without an explicit length
scale, it is chosen by log-marginal-likelihood. See the GPR section above.

Typical GPR-only config:

```
LAMBDAONLY
LAUE NORMGPR 0.780 6.300
LAUE NORMLAMREF 1.0
OUTPUT UNMERGED
```

### Typical Laue configuration

```
LAUE NORMCHEBYSHEV 6 0.780 1.500
LAUE NORMCHEBYSHEV 6 1.500 6.300
LAUE NORMLAMREF 1.0
```

Two ranges handle the characteristic intensity discontinuity at the bandpass boundary.

### Other recommended keywords for Laue data

```
SCALES CONSTANT BFACTOR OFF   # one scale per dataset, no B-factor
ANOMALOUS OFF                  # Friedel pairs merged
SDCORRECTION NOREFINE          # no SD correction refinement
```

### `LAMBDAONLY` — wavelength-normalization-only mode

```
LAMBDAONLY
LAUE NORMCHEBYSHEV 9 2.8 3.8
LAUE NORMLAMREF 3.2
OUTPUT UNMERGED
```

`LAMBDAONLY` runs **only** the wavelength normalization and nothing else — no
primary/B-factor/secondary/detector scaling, and no outlier rejection at any
stage. `SCALEUSED` in the unmerged output is then purely the wavelength scale
`ws(λ)` (see the note below on what that column means); the full observation
count is preserved (zero rejections).

Flow when active:

1. `SetOnlyLambda()` sets `initialScale = roughScale = mainScale = false` and
   `onlyLambda = true` (so `OnlyMerge()` returns false — a wavelength scale
   model is still built).
2. Outlier policy forced to `SetNoreject()` for both scaling and merging.
3. SD-correction refinement turned off (unless `SDCORRECTION REFINE` explicit).
4. The Laue pre-normalisation pass runs with the **full** `REFINE CYCLES` count
   (not the fixed 5 used in normal mode), then scales are applied (`ws` only,
   `ps=bs=ss=ds=1`) and the unmerged MTZ is written.

**No `SCALES CONSTANT` needed:** if primary scaling reports insufficient
information but LAUE wavelength ranges are present, `ScaleModel::init` falls
back to constant primary scaling automatically so only the Chebyshev
coefficients refine. (Supplying `SCALES CONSTANT` explicitly also works.)

**Keyword-name caution:** `keyIs()` compares only the first 4 characters, so the
keyword must **not** collide with `ONLYMERGE` — both `ONLY…` forms map to the
same 4-char key `ONLY` and dispatch to whichever class is inherited first in
`InputAll`. This is why the keyword is `LAMBDAONLY` (`LAMB`), not `ONLYLAMBDA`.

### The unmerged output file

`writeunmerged.cpp` writes the **scaled** intensities, not the raw ones:

```cpp
data[ic++] = this_obs.kI();      // kI() = I_/gscale   -> the SCALED intensity
float g = this_obs.Gscale();
if (g != 0.0) g = 1.0f/g;
data[ic++] = g;                  // SCALEUSED = 1/gscale
```

so `I_out = I_in × SCALEUSED`. `SCALEUSED` records the factor that **was
applied**, not one still to apply — dividing by it a second time before
computing merging statistics roughly doubles R-merge.

`LAMBDA` is written back out whenever it was read in (`data_flags::is_lambda`),
as an optional column immediately after `TIME` in the column list and in both
the raw and summed-partials data paths. Without it, scaled Laue data could not
be re-scaled or analysed against wavelength without going back to the input file
and matching on `(h, k, l, M/ISYM, batch)`. The raw path's `ic == NumCol`
assertion guards the column/data ordering if further columns are added.

### LAMBDA MTZ column

aimless reads a per-reflection wavelength column (type R) as an optional column.
Three labels are accepted, tried **in this order of preference** and each matched
**without regard to case**:

| order | label |
|-------|-------|
| 1 | `LAMBDA` |
| 2 | `LAM` |
| 3 | `WAVELENGTH` |

The first one present in the file is used; the others are ignored even if also
present. So a file carrying both `Lambda` and `wavelength` uses `Lambda`.

If none is present, each reflection inherits its batch wavelength via
`Batch::Wavelength()`. When one is present, `obs.lambda()` returns the
per-reflection wavelength used throughout the scale model.

The chosen column is reported in the log:

```
Per-reflection wavelength taken from column Lambda
```

**Implementation** (`mtz_unmerge_io.cpp`): the logical column name stays
`LAMBDA` throughout the program — only the file-to-program mapping in
`MtzUnmrgFile::get_col_lookup` is affected. `CMtz::MtzColLookup` is still tried
first (exact match); only if that fails does the alias list
`WavelengthColumnAliases` get scanned via the local helper `ColLookupNoCase`,
which walks every column of every dataset comparing `StringUtil::ToUpper`
forms. On a hit, `CNL.label` is set to the actual file label so the log and
`column_labels::Label("LAMBDA")` report what was really used. No other column
gets alias treatment.

---

## Radiation type (`PROBE NEUTRON | XRAY`)

lawless exists because LSCALE did, and LSCALE was for neutron Laue as much as
X-ray. Its entire neutron provision was one switch, `POL_Type neutron`, whose
whole effect (`Daresbury_laue/doc/lsm.txt` §2.3) is `fP = 1.0`.

```
PROBE  NEUTRON | XRAY            (default XRAY)
PROBE  NEUTRON TOF | QUASILAUE   (default TOF)
```

The second axis is the **instrument**, not the radiation: time-of-flight
separates diffraction orders in time, so there is nothing to deconvolute; a
reactor quasi-Laue instrument does not, and would need the harmonic machinery
LAUENORM has and this program does not.

`PROBE` changes **defaults, guards and reporting only**. It never applies a
correction that cannot be seen in the log or overridden by another keyword,
everything it switches stays independently settable, and the log prints a block
naming every default it changed. With `PROBE XRAY` (the default) the output is
byte-identical to the same build without the keyword.

### What `PROBE NEUTRON` changes

| | Why |
|---|---|
| `ANOMALOUS OFF`, and not switched on automatically | Nuclear scattering has no anomalous signal outside a short list of isotopes (¹¹³Cd, Sm, Gd, ¹⁰B, ⁶Li, In). Without this, aimless turns anomalous on by itself whenever it thinks it sees one. `ANOMALOUS ON` overrides |
| polarization factor set to zero | `fP = 1` for neutrons. (`UpdatePolarizationCorrections` is not called from aimless, so this is defensive) |
| `SCALES BATCH` by default | Laue exposures are stationary, so φ is not a scaling variable and rotation smoothing has nothing to smooth over. The hook is `SCALES::setDefaultBatchMode()`, which switches any specification that did not choose a mode for itself; an explicit `BATCH`, `ROTATION`, `SPACING`, `BROTATION` or `CONSTANT` sets `ScaleSpecification::scalemodegiven` and always wins. This also turns a pre-existing abort into a working default: `SCALES BFACTOR ON` with no mode falls back to rotation smoothing over a φ range that does not exist, and the refinement dies in `samplemean()` |
| radiation-damage wording becomes stability wording | Neutrons do not damage the crystal. The relative B-factor is *not* a dose correction; it is a relative resolution-dependent scale absorbing crystal slippage, centring and illuminated-volume drift |
| outlier rejections reported against 2θ | `REJECT` tests against the weighted mean, and for Laue data the weights vary systematically with λ and d, ie with 2θ. On CuZnSOD the rejection rate runs from 13.5 % in the 15–30° bin to 0.4 % at 90–105°, and switching rejection off moves the measured high-angle intensity deficit from −9.2 % to −5.1 % |

### Warnings issued after HKLIN is read

- no `LAMBDA` column, so every observation takes its wavelength from the batch
  header — fatal to wavelength normalisation, though not an error here
- batches whose `LDTYPE` is not 3 (Laue)
- `QUASILAUE`, because harmonic deconvolution is not implemented

### Regression

Adding `PROBE NEUTRON` to the established configuration for CuZnSOD, dMPro-KB5
and hCAII-Cu gives byte-identical merged and unmerged MTZs, because those decks
already set `SCALES BATCH`/`CONSTANT` and `ANOMALOUS OFF` explicitly. Dropping
those explicit keywords and relying on the `PROBE NEUTRON` defaults reproduces
the same files again on all three. With `PROBE XRAY` (the default) every output
— MTZs, log, `SCALES`, `LAMBDANORM`, `ROGUES`, `NORMPLOT`, `CORRELPLOT`,
`ANOMPLOT` — is byte-identical to the build without the keyword, once the MTZ
run timestamps are masked.

### Implementation

`Probe` (`probe.hh/.cpp`) is a **static** class, like `SelectI` in
`hkl_unmerge.hh`, so that reporting code deep in `printing.cpp` and
`radiationdamageanalysis.cpp` can ask what the radiation is without threading a
control object through every signature. It is set once in `aimless.cpp`
immediately after the `ANOMALOUS` controls, and every guard is gated on
`Probe::IsNeutron()`.

### Why an angle term is the next thing (`tools/`)

The four scripts in `tools/` are the measurement `PROBE` was built on. For
every unique reflection with ≥3 observations they compare each observation with
the leave-one-out mean of its symmetry mates; resolution is fixed within a
reflection, so a purely resolution-dependent error cancels and what survives is
a genuine systematic.

On CuZnSOD, dMPro-KB5 and hCAII-Cu the residual bias is a reproducible function
of **scattering angle**, about 9 % peak to peak, ~10σ against the null model in
`null_test.py`. A cubic in 2θ explains 53 %, 53 % and 24 % of the weighted bias
variance where a cubic in λ explains 31 %, 13 % and 33 %, and adding 2θ on top
of λ takes all three to 63–74 %. `SECONDARY 6` with real UB matrices moves the
high-angle deficit only from −9.2 % to −7.6 %: the spherical harmonics are in
the **crystal** frame, and a stationary exposure maps a lab-frame 2θ effect onto
a different part of that surface in every pack.

A post-hoc `g(2θ)` applied to already-scaled data makes R-merge *worse* by
0.5–1 % at every Chebyshev degree tried, because the batch scales, relative
B-factors and λ curve were fitted with the systematic present and have partly
compensated for it. The term has to be refined **jointly, inside the scale
model**, and judged on χ²/CC½/maps rather than R-merge.

Run `null_test.py` before believing any of the others.

---

## Key files

| File | Role |
|------|------|
| `scaletypes.hh/.cpp` | `WavelengthChebyshevScale` (log-Chebyshev); `WavelengthGPRScale` (Eigen GP fit, fixed lookup) |
| `scalemodel.hh/.cpp` | Integrates wavelength scale; wavelength-only mode; parameter management; GPR fit driver + fixed-correction application |
| `aimless.cpp` | Laue Chebyshev pre-pass + GP pre-pass (leave-one-out log-ratio sampling) before `FC.roughScale`; `LAMBDAONLY` wiring (apply ws-only scales, force no-reject) |
| `keywords_aimless.hh/.cpp` | `LAUE` keyword parser; NORMCHEBYSHEV/NORMLAMREF/NORMGPR* disambiguation; `LAMBDAONLY` keyword; `PROBE` parser and `SCALES::setDefaultBatchMode()` |
| `probe.hh/.cpp` | `Probe` — static radiation type and instrument class, from `PROBE NEUTRON\|XRAY [TOF\|QUASILAUE]` |
| `reject.hh/.cpp` | `RejectAngleStats` — outlier rejections counted against 2θ, active under `PROBE NEUTRON` |
| `printing.cpp`, `radiationdamageanalysis.cpp` | Radiation-damage wording becomes stability wording under `PROBE NEUTRON` |
| `tools/` | Leave-one-out bias diagnostics: `twotheta_probe.py`, `separability.py`, `null_test.py`, `twotheta_sim2.py` |
| `globalcontrols_aimless.hh` | `FlowControl::SetOnlyLambda()`/`OnlyLambda()`; `OnlyMerge()` excludes onlyLambda |
| `InputAll_aimless.hh` | Inherits `LAUE`, `LAMBDAONLY` into `InputAll` |
| `hkl_unmerge.hh/.cpp` | `lambda_` on `observation_part` and `observation`; `store_part` passes lambda |
| `mtz_unmerge_io.cpp` | Reads wavelength column (LAMBDA/LAM/WAVELENGTH, case-insensitive); batch-wavelength fallback; substitutes the batch number for φ if there is no ROT column |
| `writeunmerged.cpp` | Unmerged output; writes scaled I with `SCALEUSED = 1/gscale`, and `LAMBDA` when present |
| `hkl_datatypes.hh/.cpp` | `Batch::HtoSr0` — diffraction vector, takes the per-observation wavelength for Laue; `Batch::Ldtype()` |
| `hkl_unmerge.cpp` | `CalcSecondaryBeams` and friends — pass `observation::lambda()` down to `HtoSr0` |
| `sdctypes.hh/.cpp` | `SDcorrection` — the `(lambda/lambda_ref)^SdLam` factor, applied in `Correct(observation&, Iav)` |
| `sdmodel.hh/.cpp` | `SDmodel::FitLambdaSDcorrection` — per-run wavelength exponent from binned chi^2 |
| `columnlabels.hh/.cpp` | `col_lambda` column index; `is_lambda` flag in `DataFlags` |
| `openinputfile.cpp` | Registers LAMBDA as optional MTZ column |
| `CMakeLists.txt` | Build configuration with source-tree headers and CCP4-9 dylibs |
| `CMakePresets.json` | All library paths for debug/release builds |
| `.vscode/` | VSCode tasks, settings, IntelliSense config |
| `test_data/run_aimless.sh` | Test script for `combined_sorted.mtz` |

---

## Notes

- All existing monochromatic functionality is unchanged; `ws = 1.0` when `LAUE` keyword is absent.
- `NormaliseParameters()` returns early when `nprimaryscale == 0` (e.g. `SCALES CONSTANT`) to avoid an assertion on `scalenormbatch`.
- `ScaleModel::init` builds the wavelength scale even if primary `setup()` reports `status < 0`, by falling back to constant primary scaling (enables `LAMBDAONLY` without explicit `SCALES CONSTANT`).
- `FormatSave`/`Restore` on `WavelengthChebyshevScale` handle persistence in aimless checkpoint files (note: save format stores log-coefficients `a_k`).

---

## To do

- **A scattering-angle scale term (`SCALES ... TWOTHETA <n>`)** — the next
  phase, and the one the measurement in `tools/` argues for. A Chebyshev in
  `cos 2θ` per run, entering the product in `ScaleModel::ScaleFactor`
  (`scalemodel.cpp`) as one more factor; `2θ = 2 asin(λ/2d)` and both are
  already per-observation, so no new MTZ column is needed. Offer the
  one-parameter LSM obliquity form `exp(u/cos 2θ)` alongside it for comparison,
  though the measured residual is humped (−1 % at 25°, +5 % at 45°, −5 % past
  75°) and both LSM forms are monotone. Risk: correlation with the λ curve and
  the relative B-factor, all three of which reach 2θ through *d* — normalise
  `⟨ln g⟩ = 0` over the observed range and tie toward flat, as `SECONDARY` does.

- **`SECONDARY` on by default under `PROBE NEUTRON` when UB is present** — for a
  1 mm crystal in D₂O the attenuation is μ ≈ 0.195 mm⁻¹ hydrogenous, 0.071 mm⁻¹
  perdeuterated, dominated by the incoherent cross section of hydrogen
  (80.3 barn against 0.33 barn for absorption). That is a 13 % (5 %)
  transmission spread across orientation, which is real and is `SECONDARY`'s
  job. It needs a guard that refuses to fit 48 coefficients when the batch
  headers carry an identity UMAT. Note the *wavelength dependence* of that
  absorption is negligible — the 1/v term is ~1 % of μ and transmission changes
  0.2 % across a 2.1–3.9 Å band — so a λ-dependent secondary surface is not
  worth building.

- **Update wavelength normalisation during scaling refinement** — investigate whether continuing to update the wavelength normalisation throughout determination of the scaling parameters (including the pre-pass) improves overall statistics, versus fixing `ws(λ)` after the pre-pass and holding it constant during joint refinement. (Applies to both the Chebyshev pre-pass and the new GPR pre-pass, which currently fixes `ws(λ)` after fitting.)

- **GPR refinement / hyperparameters** — the GPR (`LAUE NORMGPR`) is implemented as a fixed non-parametric pre-pass (see "Gaussian-process wavelength normalisation"). `σ_f` and a noise-inflation factor are now optimised jointly with the length scale, and merging statistics have been compared against the Chebyshev model on `ca_thio` (they agree). Possible follow-ups: full marginal-likelihood *gradient* optimisation instead of the 3-D grid search; a proper CC½/half-dataset comparison; the held-out χ²/n is still ~2, suggesting the bin variances remain slightly optimistic (samples from the same reflection are correlated across bins) — a per-reflection random effect would tighten this.

- **Fit SDLAMBDA after rebuilding, on more datasets** — the wavelength term in
  the SD correction (`SDCORRECTION SDLAMBDA`) is self-limiting on the two
  datasets it has been run on (one with a wavelength residual, one without),
  but both are MaNDi. A non-monotone residual would defeat the single power
  law; see the SDLAMBDA section for what to check first.

- **Per-run GPR on a dataset that actually drifts** — `NORMGPRPERRUN` is
  implemented and regression-clean, but the only dataset it has been run on
  (CuZnSOD, 10 packs) has no per-run spectral drift: per-run `σ_f` comes out at
  0.002–0.018 in log space against 1.66 for the global fit. It has therefore
  never been exercised against a real non-flat residual.
