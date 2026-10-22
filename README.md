# lawless: aimless with wavelength normalization

lawless is a modification of CCP4 [aimless](https://www.ccp4.ac.uk/html/aimless.html) [1]. It adds tools relevant to scaling of Laue data, including:

- wavelength normalisation `w(λ)` in the scale model, as a Chebyshev
  polynomial (`LAUE NORMCHEBYSHEV`) or a Gaussian process (`LAUE NORMGPR`);
- Lorentz correction on input, `I' = I·sin²θ/λᵖ` (`LORENTZ [<p>]`, off by
  default; `LORENTZ` alone gives `p = 0`, `LORENTZ 4` the time-of-flight form);
- a wavelength term in the SD correction (`SDCORRECTION ... SDLAMBDA`);
- a **radiation-type switch** (`PROBE NEUTRON | XRAY`) that sets neutron
  defaults (e.g. no anomalous signal checks).

If these keywords are not used, it behaves the same as aimless. It is intended as a replacement for **LSCALE/LAUENORM** from the Daresbury Laue suite, and SCALA, which is the predecessor to aimless.

### Building

lawless needs a C++11 compiler, CMake 3.21 or newer, and these libraries:

| dependency | used for | Ubuntu package | macOS source |
|---|---|---|---|
| clipper 2.1 | cell, symmetry, MTZ classes | `libclipper-dev` | CCP4 (library), headers separately |
| libccp4 8.0 | MTZ file I/O | `libccp4-dev` | CCP4 (library), headers separately |
| mmdb2 2.0 | needed by clipper | `libmmdb2-dev` | CCP4 (library), headers separately |
| cctbx | space groups | `libcctbx-dev` | conda `cctbx-base` |
| boost (headers) | | `libboost-dev` | Homebrew `boost` |
| Eigen 3.4 (headers) | GP wavelength normalisation | `libeigen3-dev` | CCP4 |
| TNT/JAMA (headers) | eigenvalues | `libtnt-dev`, `libjama-dev` | CCP4 |

Nothing from CCP4 is needed at run time except the libraries themselves. The
warning `CLIBD not set` that appears when run outside a CCP4 environment can be
ignored.

#### Linux (Ubuntu 24.04)

Everything comes from the distribution, and no CCP4 installation is needed:

```bash
sudo apt install build-essential cmake libboost-dev libeigen3-dev \
    libccp4-dev libmmdb2-dev libclipper-dev libcctbx-dev libtnt-dev libjama-dev
cmake --preset release
cmake --build --preset release              # -> build/lawless
sudo cmake --install build                  # -> /usr/local/bin/lawless
```

#### macOS

A CCP4 installation supplies the clipper, libccp4 and mmdb2 libraries, Eigen
and TNT/JAMA. It does not supply cctbx as a library, or the clipper, libccp4
and mmdb2 headers, so these come from elsewhere:

1. [CCP4](https://www.ccp4.ac.uk/download/) 9. The newest
   `/Applications/ccp4-*` is found automatically.
2. CMake and boost: `brew install cmake boost`.
3. cctbx, in a conda environment:
   `conda create -n cctbx -c conda-forge cctbx-base`. The binary links to
   `libcctbx` in this environment, so keep it.
4. The clipper, libccp4 and mmdb2 headers, by either route:
   - from the brewsci/bio tap:
     `brew install brewsci/bio/libccp4 brewsci/bio/mmdb2 brewsci/bio/clipper4coot`;
     these are found in `/opt/homebrew/include` with no further setting; or
   - from the CCP4 source tarballs at <https://ftp.ccp4.ac.uk/opensource/>:
     `libccp4-8.0.0.tar.gz`, `mmdb2-2.0.22.tar.gz` and
     `clipper-2.1.20201109.tar.gz`. Unpack all three into one directory and
     set `CCP4_SRC` to it. Nothing in them needs to be built.

```bash
export CCTBX_DIR=~/miniforge3/envs/cctbx    # wherever the env is
export CCP4_SRC=~/ccp4-headers              # tarball route only
cmake --preset release
cmake --build --preset release              # -> build/lawless
cmake --install build --prefix ~/.local     # -> ~/.local/bin/lawless
```

#### Telling CMake where things are

| variable | default | |
|---|---|---|
| `CCP4_DIR` | `$CCP4`, then the newest `/Applications/ccp4-*` | CCP4 installation |
| `CCP4_SRC` | none | directory holding the clipper, libccp4 and mmdb2 sources, unpacked tarballs or a CCP4 source checkout |
| `CCTBX_DIR` | `$CONDA_PREFIX` | conda environment with `cctbx-base` |

Each can be set as an environment variable or as `-D<name>=...` on the
`cmake` command line. Exporting them (in `~/.zshrc`, say) means they survive a
deleted `CMakeCache.txt`. CMake prints every path it used, and stops with a
list of whatever it could not find.

The `default` preset is a Debug build in the same `build/` directory. Every
build also makes `build/aimless` a symlink to `lawless`; `cmake --install`
installs `lawless` only.

## How it differs from aimless

| | aimless | lawless |
|---|---|---|
| wavelength | one per batch (`ALAMBD` in the batch header) | one per observation (`LAMBDA` / `LAM` / `WAVELENGTH` column) |
| scale model | `primary × B-factor × secondary × detector` | `… × wavelength normalisation` |
| error model | `SdFac √(σ² + SdB·I + (SdAdd·I)²)` | optional wavelength term on top |
| radiation | X-ray assumed throughout | `PROBE NEUTRON` changes defaults, guards and reporting |
| Lorentz | assumed already applied | optional `LORENTZ` correction on input, for data where it is not |
| output | unmerged MTZ | unmerged MTZ now carries `LAMBDA`; adds `LAMBDANORM.plt` |

## Keywords unique to lawless

| keyword | what it does |
|---|---|
| `LAUE NORMCHEBYSHEV <degree> <lam_min> <lam_max>` | wavelength normalisation as a refined Chebyshev polynomial, up to 5 ranges |
| `LAUE NORMGPR <lam_min> <lam_max>` | wavelength normalisation as a Gaussian process, fitted once and held fixed |
| `LAUE NORMGPRLENGTH`, `NORMGPRBINS`, `NORMGPRMATERN`, `NORMGPRPERRUN`, `NORMGPRWEIGHT` | options for the GP fit |
| `LAUE NORMLAMREF <lambda_ref>` | wavelength at which `w(λ) = 1` |
| `LAMBDAONLY` | apply the wavelength normalisation and nothing else |
| `SDCORRECTION ... SDLAMBDA [<scale>]` | wavelength term `(λ/λ_ref)^SdLam` in the SD correction |
| `PROBE NEUTRON \| XRAY [TOF \| QUASILAUE]` | neutron defaults and reporting |
| `LORENTZ [<p>]` | Lorentz correction `sin²θ/λᵖ` on input |
| `PLOT GNUPLOT \| NOGNUPLOT` | gnuplot copies of the plot files |

### Requirements for the input file

The mtz input requirements for lawless are the same as aimless: the input file must be an unmerged batch mtz sorted on H K L M/ISYM BATCH. In addition, **a per-observation wavelength column is mandatory whenever the `LAUE` keyword is given.** It may be labelled `LAMBDA`, `LAM` or `WAVELENGTH` (case insensitive); lawless will search for those columns in order. 

The most convenient method for doing this is by running unmerged data through [POINTLESS](https://www.ccp4.ac.uk/html/pointless.html); however, the current release version in CCP4 v9 does not output wavelength columns (nor any column in the input that it doesn't understand; see 'optional columns' [here](https://www.ccp4.ac.uk/html/pointless.html#input_and_output_files)). This will be fixed in the next major CCP4 release. For now, unmerged batch mtz files can be easily sorted with gemmi:

```python
import gemmi

mtz = gemmi.read_mtz_file("unmerged.mtz")
mtz.switch_to_original_hkl()   # undo any existing reduction, using M/ISYM
mtz.switch_to_asu_hkl()        # reduce to the CCP4 asymmetric unit, setting M/ISYM
mtz.sort(5)                    # sort on the first five columns: H K L M/ISYM BATCH
mtz.write_to_file("sorted.mtz")
```

This keeps every column, including the wavelength, and uses the space group
already in the file; it does no symmetry determination. The file must have
`H K L M/ISYM BATCH` as its first five columns. If it has no `M/ISYM` column,
add one set to 1 (original indices) before running this. On the test data the
result is identical, row for row, to POINTLESS run with `SPACEGROUP HKLIN`.

**Batch headers** should carry:

| field | needed for |
|---|---|
| `LDTYPE = 3` | marks the data as Laue | 
| `UMAT` (orientation) | `SCALES SECONDARY` / `ABSORPTION` |

### Wavelength normalisation

Two parameterisations are offered:

**Chebyshev** (`LAUE NORMCHEBYSHEV <degree> <lam_min> <lam_max>`). Coefficients represent `h(λ) = ln f(λ)` as a Chebyshev polynomial in `z(λ)`, mapped from `[λ_min, λ_max]` to `[−1, 1]`:

```
h(λ)  = Σ_{k=0..N} a_k T_k(z)
w(λ)  = exp( h(λ) − h(λ_ref) )
```

The log parameterisation guarantees `w > 0` for any coefficients, so the refinement cannot step into a negative scale. The coefficients are refined jointly with the rest of the scale model. Up to five non-overlapping wavelength ranges are supported, each with its own degree; reflections outside all ranges get `w = 1`.

**Gaussian process regression** (`LAUE NORMGPR`). A non-parametric alternative, fitted in a pre-pass from leave-one-out log-ratios and then held **fixed** — it contributes no refinable parameters. The length scale, signal variance and a noise inflation factor are chosen by marginal likelihood unless the length scale is set explicitly.

**Which to use.** GPR is the default: it follows structure a low-order polynomial cannot, and because it is fixed after its pre-pass it cannot trade against the other scale terms. Chebyshev normalisation is preferable when you want the wavelength term refined together with everything else, or when the band is narrow and smooth (most quasi-laue neutron sources are amenable to Chebyshev refinement). NORMGPR and NORMCHEBYSHEV are mutually exclusive.

Normalisation results are reported in the log as a table and an ASCII plot, in the XML, and as a self-contained gnuplot script in the file `LAMBDANORM.plt`.

#### Chebyshev keywords

```
LAUE NORMCHEBYSHEV <degree> <lam_min> <lam_max>     # add a wavelength range (up to 5 ranges)
LAUE NORMLAMREF <lambda_ref>                        # reference wavelength
```

#### GPR keywords

```
LAUE NORMGPR <lam_min> <lam_max>      # enable GP normalisation over this range
LAUE NORMGPRLENGTH <lengthscale>      # optional: fix GP length scale (Å); default auto
LAUE NORMGPRBINS <nbins>              # optional: number of training bins; default 50
LAUE NORMGPRMATERN                    # optional: Matérn-3/2 kernel; default squared-exp
LAUE NORMGPRPERRUN                    # optional: residual w(λ) per run, on top of the global curve
LAUE NORMGPRWEIGHT EXPECTED|MEASURED  # optional: sample weights; default EXPECTED
LAUE NORMLAMREF <lambda_ref>          # shared with Chebyshev
```

### Wavelength normalisation only

```
LAMBDAONLY
```

Fits and applies the wavelength normalisation and nothing else.
The wavelength refinement runs for the full `REFINE CYCLES` for 
Chebyshev normalization. Output is unmerged mtz. 

### Wavelength term in the SD correction

```
SDCORRECTION ... SDLAMBDA [<scale>]      σ' = σ'_aimless × (λ/λ_ref)^SdLam
```

Adds a wavelength dependence to the error model, for Laue data whose scatter
varies with λ in a way the intensity-dependent terms `SdFac`, `SdB` and
`SdAdd` cannot represent. One exponent `SdLam` is fitted per run, after outlier
rejection, within resolution bins and from scaled intensities. `λ_ref` is
`LAUE NORMLAMREF`. The optional `<scale>` multiplies the fitted exponent
(default 1).

It corrects **variance**, not bias: if mean `I/<I>` trends with λ, the
wavelength normalisation is at fault and should be fixed instead. It should
fit near-zero exponents on data with no wavelength-dependent scatter. Merging
statistics barely change; what improves is the σ of the merged data.

### Radiation type

`PROBE NEUTRON` changes some defaults, guards and reporting.

```
PROBE  NEUTRON | XRAY            (default XRAY)
PROBE  NEUTRON TOF | QUASILAUE   (default TOF)
```

`PROBE NEUTRON` sets:

| setting | why | to override |
|---|---|---|
| `ANOMALOUS OFF` | nuclear scattering has no anomalous signal outside a few resonant isotopes (¹¹³Cd, ¹⁵⁷Gd, ¹⁴⁹Sm, ¹⁵¹Eu, ¹⁰B, ⁶Li, ¹¹³In) | `ANOMALOUS ON` |
| `SCALES BATCH` | Laue exposures are stationary, so rotation angle is not a scaling variable. | any explicit mode: `SCALES ROTATION`, `SPACING`, `BROTATION` or `CONSTANT` |

and changes the reporting.

- Radiation damage wording becomes stability wording. The relative B-factor is
  not a dose correction for neutrons; it absorbs crystal slippage, centring and
  illuminated-volume drift.
- The "significant anomalous signal" warning is reworded rather than
  suppressed. In neutron data a Friedel difference is almost never anomalous
  scattering, and is better read as a sign that the wavelength normalisation
  is wrong.
- Outlier rejections are also reported against 2θ.

`PROBE NEUTRON QUASILAUE` (a reactor instrument such as LADI-III) also warns
that harmonic deconvolution is not implemented. The log lists every default
that `PROBE` changed.

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
default). This is for certain routines that do not perform Lorentz corrections 
after integration.

### Plot files

```
PLOT GNUPLOT | NOGNUPLOT      (default GNUPLOT)
```

Alongside the xmgrace/loggraph plot files inherited from aimless (`NORMPLOT`,
`ANOMPLOT`, `CORRELPLOT`, `ROGUEPLOT`), lawless writes self-contained gnuplot
(≥ 5.0) scripts `NORMPLOT.plt`, `ANOMPLOT.plt`, `CORRELPLOT.plt` and
`ROGUEPLOT.plt`, plus `LAMBDANORM.plt` for the wavelength normalisation. Open
any of them with `gnuplot -p <file>`. `NOGNUPLOT` turns them all off;
`PLOT NOXMGR` still controls the xmgrace files independently.

### References
[1] P.R. Evans and G. N. Murshudov, "How good are my data and what is the resolution?", *Acta Cryst.* (2013) **D69**, 1204-1214. doi: [10.1107/S0907444913000061](https://doi.org/10.1107/S0907444913000061)