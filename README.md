# lawless: aimless with wavelength normalization

lawless is a modification of CCP4 aimless. It adds tools relevant to scaling of Laue data, including:

- **wavelength normalisation** `w(λ)` in the scale model, as a Chebyshev
  polynomial (`LAUE NORMCHEBYSHEV`) or a Gaussian process (`LAUE NORMGPR`);
- **Lorentz correction** on input, `I' = I·sin²θ/λᵖ` (`LORENTZ [<p>]`, off by
  default; `LORENTZ` alone gives `p = 0`, `LORENTZ 4` the time-of-flight form);
- a **wavelength term in the SD correction** (`SDCORRECTION ... SDLAMBDA`);
- a **radiation-type switch** (`PROBE NEUTRON | XRAY`) that sets neutron
  defaults and guards.

If these keywords are not used, it behaves the same as aimless. It is intended as a replacement for **LSCALE/LAUENORM** from the Daresbury Laue suite, and SCALA, which is the predecessor to aimless. It is tailored toward neutron Laue experiments

### Building
```bash
cmake --preset default        # configure (uses CMakePresets.json)
make -C build -j4             # compile -> build/lawless
```

CMakePresets.json in the project root defines all library paths.

### Where it differs from aimless

| | aimless | lawless |
|---|---|---|
| wavelength | one per batch (`ALAMBD` in the batch header) | one per observation (`LAMBDA` / `LAM` / `WAVELENGTH` column) |
| scale model | `primary × B-factor × secondary × detector` | `… × wavelength normalisation` |
| error model | `SdFac √(σ² + SdB·I + (SdAdd·I)²)` | optional wavelength term on top |
| radiation | X-ray assumed throughout | `PROBE NEUTRON` changes defaults, guards and reporting |
| Lorentz | assumed already applied | optional `LORENTZ` correction on input, for data where it is not |
| output | unmerged MTZ | unmerged MTZ now carries `LAMBDA`; adds `LAMBDANORM.plt` |

### Changed requirements for the input file

The mtz input requirements for lawless are the same as aimless: the input file must be an unmerged batch mtz sorted on H K L M/ISYM BATCH. In addition, **a per-observation wavelength column is mandatory whenever the `LAUE` keyword is given.** It may be labelled `LAMBDA`, `LAM` or `WAVELENGTH` (case insensitive); lawless will search for those columns in order. 

The most convenient method for doing this is by running unmerged data through **POINTLESS**, however, the current version does not output wavelength columns. This will be fixed in the next major CCP4 release. 

**Batch headers** should carry:

| field | needed for |
|---|---|
| `LDTYPE = 3` | marks the data as Laue | 
| `UMAT` (orientation) | `SCALES SECONDARY` / `ABSORPTION` |
| `ROT` column (optional) | per-batch φ |

### The scale model

For each observation, the scale applied is a product:

```
g = g_primary(batch or φ) × g_B(batch, 1/d²) × g_secondary(ŝ₂) × g_detector(x,y) × w(λ)
```

The first four are identical to aimless. `w(λ)` is the wavelength normalisation described in the next section. The corrected intensity is `I / g`. 

The Lorentz correction (`LORENTZ`) is not part of `g`. When requested it is applied to the intensities and their standard deviations as the file is read, before any scaling, and it is not refined.

### Wavelength normalisation

The measured intensity of a Laue reflection is a function of the probe wavelength, along with the wavelength dependence of the detector efficiency and of everything else in the beam. `w(λ)` absorbs all of this into a single empirical curve, normalised so that `w(λ_ref) = 1`.

Two parameterisations are offered:

**Chebyshev** (`LAUE NORMCHEBYSHEV <degree> <lam_min> <lam_max>`). Coefficients represent `h(λ) = ln f(λ)` as a Chebyshev polynomial in `z(λ)`, mapped from `[λ_min, λ_max]` to `[−1, 1]`:

```
h(λ)  = Σ_{k=0..N} a_k T_k(z)
w(λ)  = exp( h(λ) − h(λ_ref) )
```

The log parameterisation guarantees `w > 0` for any coefficients, so the refinement cannot step into a negative scale. The coefficients are refined jointly with the rest of the scale model. Up to five non-overlapping wavelength ranges are supported, each with its own degree; reflections outside all ranges get `w = 1`.

**Gaussian process regression** (`LAUE NORMGPR`). A non-parametric alternative, fitted in a pre-pass from leave-one-out log-ratios and then held **fixed** — it contributes no refinable parameters. The length scale, signal variance and a noise inflation factor are chosen by marginal likelihood unless the length scale is set explicitly.

**Which to use.** The GPR is the better default: it follows structure a low-order polynomial cannot, and because it is fixed after its pre-pass it cannot trade against the other scale terms. Chebyshev normalisation is preferable when you want the wavelength term refined together with everything else, or when the band is narrow and smooth (most quasi-laue neutron sources are amenable to Chebyshev refinement). NORMGPR and NORMCHEBYSHEV are mutually exclusive.

Normalisation results are reported in the log as a table and an ASCII plot, in the XML, and as a self-contained gnuplot script in the file `LAMBDANORM.plt`.

**Keywords**

```
LAUE NORMCHEBYSHEV <degree> <lam_min> <lam_max>
LAUE NORMLAMREF <lambda_ref>
```

- `NORMCHEBYSHEV` — add a wavelength range (repeat for each range, up to 5)
- `NORMLAMREF` — reference wavelength; `ws(λ_ref) = 1.0` by definition

#### GPR (Gaussian-process) alternative

```
LAUE NORMGPR <lam_min> <lam_max>   # enable GP normalisation over this range
LAUE NORMGPRLENGTH <lengthscale>   # optional: fix GP length scale (Å); default auto
LAUE NORMGPRBINS <nbins>           # optional: number of training bins; default 50
LAUE NORMGPRMATERN                 # optional: Matérn-3/2 kernel; default squared-exp
LAUE NORMGPRPERRUN                 # optional: residual w(λ) per run, on top of the global curve
LAUE NORMLAMREF <lambda_ref>       # shared with Chebyshev
```

### Radiation type

`PROBE NEUTRON` changes **defaults, guards and reporting only**. It never applies a correction that cannot be seen in the log or overridden by another keyword, every term it switches stays independently settable, and the log prints a block naming exactly what it changed. With `PROBE XRAY` — the default — the output is identical to the same build without the keyword.

### Lorentz correction

```
LORENTZ [<p>]

LORENTZ          I' = I * sin^2(theta)             geometric Laue factor (p = 0)
LORENTZ 4        I' = I * sin^2(theta) / lambda^4  time-of-flight convention
LORENTZ <p>      I' = I * sin^2(theta) / lambda^p
```

Applies the Laue Lorentz correction to the intensities **on input, before
anything else is done with them**: to `I` and `SIGI`, and to the profile-fitted
pair `IPR` and `SIGIPR` when present. The optional number `p` is the power of
wavelength; it defaults to 0. Without the keyword no correction is made (the
default). Default is no Lorentz correction. This is for certain routines that
do not perform Lorentz corrections after integration.
