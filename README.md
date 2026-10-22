# lawless: aimless with wavelength normalization

lawless is a modification of CCP4 [aimless](https://www.ccp4.ac.uk/html/aimless.html) \[1\]. It adds tools relevant to scaling of Laue data, including:

- **wavelength normalisation** `w(λ)` in the scale model, as a Chebyshev
  polynomial (`LAUE NORMCHEBYSHEV`) or a Gaussian process (`LAUE NORMGPR`);
- **Lorentz correction** on input, `I' = I·sin²θ/λᵖ` (`LORENTZ [<p>]`, off by
  default; `LORENTZ` alone gives `p = 0`, `LORENTZ 4` the time-of-flight form);
- a **wavelength term in the SD correction** (`SDCORRECTION ... SDLAMBDA`);
- a **radiation-type switch** (`PROBE NEUTRON | XRAY`) that sets neutron
  defaults and guards.

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

The Xcode compiler has no OpenMP, so the macOS build runs on one thread.

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

### How it differs from aimless

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

### References
[1] P.R. Evans and G. N. Murshudov, "How good are my data and what is the resolution?", *Acta Cryst.* (2013) **D69**, 1204-1214. doi: [10.1107/S0907444913000061](https://doi.org/10.1107/S0907444913000061)