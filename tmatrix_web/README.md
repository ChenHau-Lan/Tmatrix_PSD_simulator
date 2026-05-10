# Local T-matrix PSD Explorer

This folder reorganizes the workspace into a myPSD-like web app:

- `backend/app/local_tmatrix.py` builds and caches local Fortran T-matrix scatter tables, then integrates PSDs into radar variables.
- `backend/app/main.py` serves the frontend and exposes `POST /api/compute`.
- `frontend/index.html` provides the PSD controls, plot, metrics table, and modelling assumptions.

Run it from the repository root:

```bash
python3 -m tmatrix_web.backend.app.main
```

Then open:

```text
http://127.0.0.1:7860
```

## Difference from snesbitt/myPSD

The public myPSD app uses a modern Docker/FastAPI/React stack and computes
T-matrix tables through `rustmatrix`, a Rust port of `pytmatrix`. It supports
rain, hail, and several snow habits with per-species assumptions, cached scatter
tables, and a JSON API.

This local app keeps the same separation of concerns, but the scattering engine
is the Fortran code already present in `py_Tmatrix_Mueller/fortran_tm`. Rain is
the verified path in this workspace. Other hydrometeor choices are exposed
through the legacy `dclass` switches and are marked experimental in the UI until
their refractive index, density, axis-ratio, and validation cases are checked.

Two important local fixes were required:

- The legacy wrapper used shell strings to run the executables, which breaks
  when the repository path contains a space (`radar simulator`).
  `py_Tmatrix_Mueller/scatter.py` now uses `subprocess.run` with argument
  lists, so the executables run correctly from this folder.
- The legacy Mueller executable opens hard-coded relative files (`inp`,
  `inp1`, `inp2`) for non-rain species. The wrapper now keeps those aliases in
  sync with the freshly generated `out1_tmat`, so hail/graupel/snow diagnostic
  runs no longer fail at EOF.

## Hydrometeor Validation Status

Rain is executable and previously verified in this workspace.

For rain, the local app now mirrors the myPSD DSD handling:

- The PSD curve is evaluated on `D = 0.05-19.95 mm`, matching the myPSD plot
  grid and total-number integration.
- Radar-variable integration uses the rain scatter table through `D_max =
  10 mm`, matching myPSD's rain scatter-table bound.
- The legacy `Scatterer.radius` name is treated as equivalent diameter in mm,
  because the bundled examples pass `aD` directly and the Fortran input writes
  that value as `DEQ` in cm. Passing `D/2` was the main cause of the earlier
  rain mismatch.

Frozen rain, dry hail, spongy graupel, and spongy aggregate now execute through
the legacy species switches, but they are not scientifically validated yet. The
source file `py_tmat.f90` contains documented dielectric constants for dry ice
and low-density aggregates, and this repo now includes source edits to apply
those constants by `IOPT`. However, the local `gfortran` toolchain currently
fails at link time with `ld: library not found for -lSystem`, so the bundled
executables could not be rebuilt in this session. Until they are rebuilt and
checked against reference cases, non-rain output should be treated as a
diagnostic of the plumbing rather than as trusted hydrometeor physics.
