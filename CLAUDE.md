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
| **Wavelength** | `ws` | Chebyshev normalisation (Laue only) |

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
stage. The unmerged output contains the original intensities with `SCALEUSED`
set purely from the wavelength scale `ws(λ)`; the full observation count is
preserved (zero rejections).

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

### LAMBDA MTZ column

aimless reads `LAMBDA` (type R) as an optional column. If absent, each reflection inherits its batch wavelength via `Batch::Wavelength()`. With LAMBDA present, `obs.lambda()` returns the per-reflection wavelength used throughout the scale model.

---

## Key files

| File | Role |
|------|------|
| `scaletypes.hh/.cpp` | `WavelengthChebyshevScale` class; log-Chebyshev algorithm |
| `scalemodel.hh/.cpp` | Integrates wavelength scale; wavelength-only mode; parameter management |
| `aimless.cpp` | Laue pre-normalisation pass (before `FC.roughScale`); `LAMBDAONLY` wiring (apply ws-only scales, force no-reject) |
| `keywords_aimless.hh/.cpp` | `LAUE` keyword parser; NORMCHEBYSHEV/NORMLAMREF disambiguation; `LAMBDAONLY` keyword |
| `globalcontrols_aimless.hh` | `FlowControl::SetOnlyLambda()`/`OnlyLambda()`; `OnlyMerge()` excludes onlyLambda |
| `InputAll_aimless.hh` | Inherits `LAUE`, `LAMBDAONLY` into `InputAll` |
| `hkl_unmerge.hh/.cpp` | `lambda_` on `observation_part` and `observation`; `store_part` passes lambda |
| `mtz_unmerge_io.cpp` | Reads LAMBDA column; batch-wavelength fallback |
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

- **GPR wavelength normalisation** — replace or supplement the Chebyshev pre-pass with a Gaussian Process smoother. Eigen3 (already in include path) provides the required Cholesky decomposition. Intended as a non-parametric pre-pass: bin reflections by λ, compute normalised bin means, fit GP with squared-exponential or Matérn-3/2 kernel, apply correction as a fixed lookup table before iterative scaling. Hyperparameter selection via marginal likelihood or heuristic length-scale. Advantage over Chebyshev: naturally adapts to variable data density and avoids polynomial runaway at sparse wavelength extremes.
