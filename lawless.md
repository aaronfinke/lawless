# LAWLESS (CCP4: Pre-release Program)

> **Draft for review.** This is the source for `lawless.html`, which will replace `aimless.html` for the Laue branch. Sections marked **[inherited]** carry over from `aimless.html` unchanged and are not reproduced here — only their heading is shown, so the shape of the final document is visible. Everything else is new or changed and is written out in full.

---

## NAME

**lawless** — scale together multiple observations of reflections from Laue data, with wavelength normalisation

## SYNOPSIS

```
lawless HKLIN foo_in.mtz HKLOUT foo_out.mtz
```

[Keyworded Input](#keyworded-input) · [References](#references)

---

## DESCRIPTION

**lawless** is the Laue branch of **aimless**. It retains the full features of **aimless**: scaling, outlier rejection, error-model correction, merging, and full statistical analysis; and adds tasks specific to Laue data:

- a **wavelength normalisation** term in the scale model, which corrects for the incident spectrum and the wavelength dependence of the detector, fitted either as a Chebyshev polynomial or as a Gaussian process;
- a **radiation-type switch**, `PROBE`, which sets defaults and guards appropriate to neutrons.

Monochromatic data are scaled exactly as by aimless: with no `LAUE` keyword the wavelength term is absent and the program behaves identically.

Lawless is intended as a replacement for **LSCALE** from the Daresbury Laue suite, and in particular for neutron Laue work.

### What is different from aimless

| | aimless | lawless |
|---|---|---|
| wavelength | one per batch (`ALAMBD` in the batch header) | one per observation (`LAMBDA` / `LAM` / `WAVELENGTH` column) |
| scale model | `primary × B-factor × secondary × detector` | `… × wavelength normalisation` |
| error model | `SdFac √(σ² + SdB·I + (SdAdd·I)²)` | optional wavelength term on top |
| radiation | X-ray assumed throughout | `PROBE NEUTRON` changes defaults, guards and reporting |
| output | unmerged MTZ | unmerged MTZ now carries `LAMBDA`; adds `LAMBDANORM` |

### Requirements for the input file

The input must be an unmerged MTZ **sorted by POINTLESS**, with two extra requirements. Note that the version of POINTLESS distributed by CCP4 (as of 09/2026) does not pass wavelength columns through to its output, nor other columns it does not itself use; a patch will be rolled out soon. 

**A per-observation wavelength column is mandatory whenever the `LAUE` keyword is given.** It may be labelled `LAMBDA`, `LAM` or `WAVELENGTH` (case insensitive). If more than one wavelength column is present, the first in that order is used and the others are ignored; for example, a file carrying both `Lambda` and `wavelength` uses `Lambda`. The label actually found is reported in the log:

```
Per-reflection wavelength taken from column Lambda
```

If the LAUE keyword is given but the input file has no wavelength column, lawless stops.

**Batch headers** should carry:

| field | needed for | if absent |
|---|---|---|
| `LDTYPE = 3` | marks the data as Laue | warning only |
| `UMAT` (orientation) | `SCALES SECONDARY` / `ABSORPTION` | the correction refines flat and does nothing — no warning is given |
| `ROT` column | per-batch φ | the batch number is substituted for φ |

**A note on POINTLESS.** It cannot reliably determine the Laue group from *unnormalised* Laue data: the wavelength dependence inflates the disagreement between symmetry mates, so the identity operation itself scores poorly and the ranking is unreliable. Give the answer explicitly:

```
LAUEGROUP P 6/m m m
CHOOSE SPACEGROUP P 65 2 2
```

### The scale model

For each observation, the scale applied is a product:

```
g = g_primary(batch or φ) × g_B(batch, 1/d²) × g_secondary(ŝ₂) × g_detector(x,y) × w(λ)
```

The first four are as described in `aimless.html`. `w(λ)` is the wavelength normalisation described in the next section. The corrected intensity is `I / g`.

### Wavelength normalisation

The measured intensity of a Laue reflection is a function of the probe wavelength, together with the wavelength dependence of the detector efficiency and of everything else in the beam. `w(λ)` absorbs all of this into a single empirical curve, normalised so that `w(λ_ref) = 1`.

Two parameterisations are offered.

**Chebyshev** (`LAUE NORMCHEBYSHEV`). Coefficients represent `h(λ) = ln f(λ)` as a Chebyshev polynomial in `z(λ)`, mapped from `[λ_min, λ_max]` to `[−1, 1]`:

```
h(λ)  = Σ_{k=0..N} a_k T_k(z)
w(λ)  = exp( h(λ) − h(λ_ref) )
```

The log parameterisation guarantees `w > 0` for any coefficients, so the refinement cannot step into a negative scale. The coefficients are refined jointly with the rest of the scale model. Up to five non-overlapping wavelength ranges are supported, each with its own degree; reflections outside all ranges get `w = 1`.

**Gaussian process regression** (`LAUE NORMGPR`). A non-parametric alternative, fitted in a pre-pass from leave-one-out log-ratios and then held **fixed** — it contributes no refinable parameters. The length scale, signal variance and a noise inflation factor are chosen by marginal likelihood unless the length scale is set explicitly.

**Which to use.** The GPR is the better default: it follows structure a low-order polynomial cannot, and because it is fixed after its pre-pass it cannot trade against the other scale terms. Chebyshev normalisation is preferable when you want the wavelength term refined together with everything else, or when the band is narrow and smooth. NORMGPR and NORMCHEBYSHEV are mutually exclusive.

Normalisation results are reported in the log as a table and an ASCII plot, in the XML, and as a self-contained gnuplot script in the file `LAMBDANORM`.

### Radiation type

`PROBE NEUTRON` changes **defaults, guards and reporting only**. It never applies a correction that cannot be seen in the log or overridden by another keyword, every term it switches stays independently settable, and the log prints a block naming exactly what it changed. With `PROBE XRAY` — the default — the output is identical to the same build without the keyword.

See the [PROBE](#probe) keyword below for the full list.

### Workflow for Laue data

1. **Integrate**, applying the Lorentz factor and any absorption correction at that stage, as usual. lawless does not apply a Lorentz factor.
2. **Sort and set symmetry with POINTLESS**, supplying `LAUEGROUP` and `CHOOSE SPACEGROUP` explicitly (see above). Check that the wavelength column survives into the sorted file — older versions of POINTLESS drop it, and lawless will then refuse to run. See [Requirements for the input file](#requirements-for-the-input-file).
3. **Scale with lawless.** A reasonable starting point for a neutron time-of-flight Laue dataset:

```
PROBE NEUTRON
RUN 1 BATCH 1 TO 11
SCALES BATCH BFACTOR ON SECONDARY 0
LAUE NORMGPR 2.10 3.90
LAUE NORMLAMREF 3.00
SDCORRECTION NOREFINE 1.0 0.0 0.0
REJECT 4
OUTPUT MERGED UNMERGED
```

4. **Check the normalisation curve** in `LAMBDANORM` before believing the statistics. A curve with structure finer than the physics can justify usually means too many GP training bins.
5. **Check χ²** by resolution and intensity, and the analysis against batch. See the notes on the error model below.

### Sections inherited from aimless.html

The following describe behaviour lawless shares with aimless and are reproduced there unchanged:

- **[inherited]** Running the program
- **[inherited]** Scaling options
- **[inherited]** Control of flow through the program
- **[inherited]** Partially recorded reflections
- **[inherited]** Scaling algorithm
- **[inherited]** Scaling to reference
- **[inherited]** Data from Denzo
- **[inherited]** Datasets

---

## KEYWORDED INPUT

Keywords new to lawless, or whose behaviour differs from aimless, are given in full. The rest are unchanged and are listed at the end.

### LAUE

```
LAUE NORMCHEBYSHEV <degree> <lam_min> <lam_max>
LAUE NORMLAMREF    <lambda_ref>
LAUE NORMGPR       <lam_min> <lam_max>
LAUE NORMGPRLENGTH <lengthscale>
LAUE NORMGPRBINS   <nbins>
LAUE NORMGPRMATERN
LAUE NORMGPRPERRUN
```

Enables Laue wavelength normalisation. Each sub-keyword may be given on its own `LAUE` line or several on one line.

| sub-keyword | meaning |
|---|---|
| `NORMCHEBYSHEV <degree> <lam_min> <lam_max>` | add a wavelength range fitted by a Chebyshev polynomial of the given degree. Repeat for up to 5 non-overlapping ranges. |
| `NORMLAMREF <lambda_ref>` | reference wavelength, at which `w = 1` by definition. Choose a wavelength near the centre of the band where the data are strong. |
| `NORMGPR <lam_min> <lam_max>` | fit the normalisation as a Gaussian process over this range instead. Values ≤ 0 take the limits from the data. |
| `NORMGPRLENGTH <lengthscale>` | fix the GP length scale, in Å. Omit (or give ≤ 0) to optimise it by marginal likelihood. |
| `NORMGPRBINS <nbins>` | number of wavelength bins in the GP training set. Default **50**. Reduce it if the log warns that the length scale has reached its lower search limit. |
| `NORMGPRMATERN` | use a Matérn-3/2 kernel instead of the default squared-exponential. Less smooth; appropriate if the true curve has kinks. |
| `NORMGPRPERRUN` | **Experimental.** After the global GP fit, fit a residual `w(λ)` for each run. It has never been exercised on a dataset with real per-run spectral drift. Use only when runs are expected to differ spectrally — on a stable source they will not, and the per-run curves come out flat. |

`NORMCHEBYSHEV` and `NORMGPR` are mutually exclusive keywords; giving both is a fatal error.

**A per-observation wavelength column in HKLIN is required** — `LAMBDA`, `LAM` or `WAVELENGTH`, case-insensitive. See [Requirements on the input file](#requirements-on-the-input-file).

### PROBE

```
PROBE NEUTRON | XRAY   [TOF | QUASILAUE]
```

Radiation type and instrument class. Default `XRAY`; for `NEUTRON` the instrument class defaults to `TOF`.

The second axis describes the **instrument**, not the radiation: time-of-flight separates diffraction orders in time, so there are no harmonics to deconvolute; a reactor quasi-Laue instrument does not separate them, and would need harmonic deconvolution, which this program does not implement.

`PROBE NEUTRON` changes the following, and prints each change it made:

| | why |
|---|---|
| `ANOMALOUS` is off, and is **not** switched on automatically | Neutron scattering lengths have no imaginary component except in a few isotopes with a nuclear resonance near thermal energies (¹¹³Cd, ¹⁵⁷Gd, ¹⁴⁹Sm, ¹⁵¹Eu, ¹⁰B, ⁶Li, ¹¹³In). Without this, aimless switches anomalous on by itself whenever it thinks it sees a signal. `ANOMALOUS ON` overrides. |
| the polarisation factor is set to zero | There is no polarisation factor for neutrons — `fP = 1`. |
| `SCALES BATCH` becomes the default | Laue exposures are stationary, so φ is not a scaling variable and rotation smoothing has nothing to smooth over. Any explicit `BATCH`, `ROTATION`, `SPACING`, `BROTATION` or `CONSTANT` overrides this. |
| radiation-damage wording becomes stability wording | Neutrons do not cause radiation damage like X-rays. The relative B-factor is not a dose correction; it is a relative resolution-dependent scale absorbing crystal slippage, centring and changes in illuminated volume. |
| outlier rejections are reported against 2θ | `REJECT` tests against a weighted mean whose weights vary with wavelength and resolution, so for Laue data the rejection rate is a function of scattering angle. The table makes that visible. |
| the "significant anomalous signal" warning is reworded | It is almost never an anomalous signal. |


### LAMBDAONLY

```
LAMBDAONLY
```

Fit and apply the wavelength normalisation only: no primary scales, no B-factors, no secondary correction, and no outlier rejection at any stage. The SD-correction refinement is also switched off unless `SDCORRECTION REFINE` is given explicitly. Useful for isolating the wavelength term, or for preparing data for a program that will do its own scaling.

### SDCORRECTION — wavelength term

The `SDCORRECTION` keyword is otherwise **[inherited]** from `aimless.html`. lawless adds one sub-keyword:

```
SDCORRECTION ... SDLAMBDA [<scale>]
```

Adds a wavelength-dependent factor to the corrected standard deviation:

```
σ' = SdFac · (λ/λ_ref)^SdLam · √(σ² + SdB·I + (SdAdd·I)²)
```

`SdLam` is **fitted from the data**, not given: lawless bins the normalised deviations by resolution and wavelength, fits the slope of log χ² against log(λ/λ_ref), and takes half of it. The optional `<scale>` multiplies the fitted exponent — `SDLAMBDA 0` disables the term, `SDLAMBDA 0.5` applies half of it. The exponent is clamped to ±2.

The fit is made **after** merge-stage outlier rejection, so it describes the data that will actually be merged.

Use it when χ² shows a trend against wavelength that the scale model cannot remove — a variance, not a bias, so no scale term can touch it. It is self-limiting: on a dataset with no wavelength trend it fits an exponent near zero and does nothing.

### SCALES

**[inherited]**, with one addition: under `PROBE NEUTRON`, `BATCH` becomes the default primary scale mode when no mode is given explicitly. See [PROBE](#probe).

### ANOMALOUS

**[inherited]**, with one addition: under `PROBE NEUTRON` the default is off and the automatic switch-on is disabled. See [PROBE](#probe).

### Unchanged keywords

The following behave exactly as documented in `aimless.html`:

**[inherited]** `RUN` · `USESDPARAMETER` · `PARTIALS` · `INITIAL` · `INTENSITIES` · `REJECT` · `ICERING` · `RESOLUTION` · `TITLE` · `ANALYSIS` · `ONLYMERGE` · `DUMP` · `RESTORE` · `REFINE` · `EXCLUDEBATCH` · `TIE` · `OUTPUT` · `UNMERGEDOUT` · `KEEP` · `LINK` · `UNLINK` · `BINS` · `NAME` · `HKLIN` · `HKLOUT` · `XMLOUT` · `HKLREF` · `LABREF` · `XYZIN` · `ROGUES` · `PLOT` · `CELL` · `BFACTOR`

---

## INPUT AND OUTPUT FILES

### Input

**[inherited]**, with one addition:

| column | | |
|---|---|---|
| `LAMBDA`, `LAM` or `WAVELENGTH` | wavelength of each observation, in Ångström. Matched without regard to case; the first present in that order is used | **required** with `LAUE`, otherwise optional |

### Output

#### Reflection files output

**[inherited]**, with one addition: the unmerged output file (`UNMERGEDOUT`) now carries a `LAMBDA` column whenever the input had a wavelength column — always written out under the logical name `LAMBDA`, whatever it was labelled on input — so the scaled data can be re-analysed against wavelength without going back to the input.

Note that intensities in the unmerged output are **scaled** — `I_out = I_in × SCALEUSED`, where `SCALEUSED = 1/g`. Do not divide by `SCALEUSED` again.

#### Other output files

**[inherited]**, with one addition:

| file | |
|---|---|
| `LAMBDANORM` | self-contained gnuplot script drawing the fitted wavelength normalisation curve. The GP fit carries a 1σ band; the Chebyshev fit does not, having no posterior covariance. Open with `gnuplot -p LAMBDANORM`. |

---

## REFERENCES

**[inherited]** — the aimless references, plus:

- Arzt, S., Campbell, J.W., Harding, M.M., Hao, Q. & Helliwell, J.R. (1999). *LSCALE — the new normalization, scaling and absorption correction program in the Daresbury Laue software suite.* J. Appl. Cryst. **32**, 554–562.
- Rasmussen, C.E. & Williams, C.K.I. (2006). *Gaussian Processes for Machine Learning.* MIT Press. — for the GP normalisation.

---

## Appendices

- **[inherited]** Appendix 1: Partially recorded reflections
- **[inherited]** Appendix 2 onwards, as in `aimless.html`

