# Changelog

User-visible changes to **lawless**, newest first. lawless is the Laue branch of
aimless; changes to aimless itself up to the fork (0.8.3) are in the release
notes in `aimless.html`.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Each version is tagged `v<version>` in git, and `version.hh` carries the same
number.

## [Unreleased]

### Fixed
- **`SCALES BATCH` no longer aborts on runs with a single batch.** Laue files
  in which every exposure (pack) has its own non-contiguous batch number, as
  DIALS writes for MaNDi, are split into one run per batch, and BATCH scaling
  stopped with "run 1 has no useful batch information for BATCH scaling". A
  single-batch run now gets one scale (and B-factor, if requested) for the
  run, as `SCALES CONSTANT` would give it, and the log says so for each run.
  This also makes the `PROBE NEUTRON` default (`SCALES BATCH`) work on such
  files.
- **`SECONDARY` and `ABSORPTION` no longer require a non-zero batch phi
  range.** The secondary-beam geometry uses only each observation's `ROT` and
  the batch orientation, so stills, which legitimately have a zero range, were
  refused for no reason. They now need a `ROT` column and valid orientation
  instead; without `ROT`, phi is taken from the batch number and the geometry
  would be wrong.

## [0.1.2] - 2026-10-22

### Changed
- **The build no longer depends on paths on the author's machine.** CMake now
  searches for every dependency. On Ubuntu 24.04 it builds from distribution
  packages alone. On macOS it uses the libraries from a CCP4 installation,
  cctbx from conda and boost from Homebrew; the clipper, libccp4 and mmdb2
  headers, which CCP4 does not install, come from Homebrew (brewsci/bio) or
  from the CCP4 source tarballs. Locations can be given with `CCP4_DIR`,
  `CCP4_SRC` and `CCTBX_DIR`. See "Building" in `README.md`.
- CMake 3.16 or newer is required (3.21 to use the presets).

### Fixed
- `cmake --install` on macOS produced a binary that could not find
  `libcctbx.dylib`.
- Builds against recent boost (1.92 tested), where `boost/timer.hpp` is
  otherwise an error.

## [0.1.1] - 2026-10-21

### Changed
- **GP wavelength normalisation now weights its training samples by the
  expected standard deviation, not the measured one** (`LAUE NORMGPRWEIGHT
  EXPECTED`, the new default). Where σ grows with I, measured weights favour
  observations that came out low by chance, so the binned ratios became
  harmonic means biased low and `w(λ)` oscillated. On simulated data this gave
  CC½ below 0.4; with expected weights it is 0.95–0.98. On real data from
  MaNDi, LADI-III and others the results are essentially unchanged.
- The binary is now named `lawless` (`build/lawless`). A `build/aimless` symlink
  is made for local compatibility.

### Added
- `LAUE NORMGPRWEIGHT EXPECTED | MEASURED`. `MEASURED` reproduces the 0.1.0
  weights exactly. The default rarely needs to be changed.
- A log warning when the GP noise inflation factor reaches the upper limit of
  its search grid.

### Documentation
- `lawless.md`: `LORENTZ` described in the overview, the scale model and the
  unmerged output; `NORMGPRWEIGHT` added; parameters that rarely need changing
  are marked as such.

## [0.1.0] - 2026-10-19

First release of lawless, forked from aimless 0.8.3. With no `LAUE` keyword,
monochromatic data are scaled exactly as by aimless.

### Added
- **Wavelength normalisation** `w(λ)` as a factor in the scale model, through
  the `LAUE` keyword:
  - `NORMCHEBYSHEV <degree> <lam_min> <lam_max>`: a Chebyshev polynomial in
    log space, up to 5 ranges, refined with the rest of the scale model after
    a wavelength-only pre-pass;
  - `NORMGPR <lam_min> <lam_max>`: a Gaussian process fitted once in a
    pre-pass and then held fixed, with `NORMGPRLENGTH`, `NORMGPRBINS` and
    `NORMGPRMATERN` options;
  - `NORMGPRPERRUN` (experimental): a residual `w(λ)` for each run on top of
    the global curve;
  - `NORMLAMREF <lambda_ref>`: the wavelength at which `w = 1`.
- `LAMBDAONLY`: fit and apply the wavelength normalisation only, with no other
  scaling and no outlier rejection.
- A per-observation wavelength column, read from `LAMBDA`, `LAM` or
  `WAVELENGTH` (any case). It is required with `LAUE`, and missing it is now a
  fatal error.
- `SDCORRECTION ... SDLAMBDA [<scale>]`: a fitted wavelength term
  `(λ/λ_ref)^SdLam` in the SD correction.
- `PROBE NEUTRON | XRAY [TOF | QUASILAUE]`: neutron defaults and guards.
  Anomalous is off and not switched on automatically, `SCALES BATCH` is the
  default, radiation-damage wording becomes stability wording, the anomalous
  signal warning is reworded, and outlier rejections are reported against 2θ.
  Every change is printed in the log.
- `LORENTZ [<p>]`: Lorentz correction `sin²θ/λᵖ` applied on input. Off unless
  asked for.
- Output of the wavelength normalisation: a table and ASCII plot in the log, a
  block in XMLOUT, and a gnuplot script `LAMBDANORM.plt`. The GP version also
  has a 1σ band and the binned training data.
- gnuplot copies (`.plt`) of `NORMPLOT`, `ANOMPLOT`, `CORRELPLOT` and
  `ROGUEPLOT`, switched by `PLOT [NO]GNUPLOT`.
- Documentation in `lawless.md`.

### Changed
- The secondary-beam direction is computed from each observation's own
  wavelength instead of the batch wavelength. This matters for `SECONDARY` and
  `ABSORPTION` on Laue data.
- The unmerged output file carries a `LAMBDA` column when the input had one.

[Unreleased]: https://github.com/aaronfinke/lawless/compare/v0.1.2...HEAD
[0.1.2]: https://github.com/aaronfinke/lawless/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/aaronfinke/lawless/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/aaronfinke/lawless/releases/tag/v0.1.0
