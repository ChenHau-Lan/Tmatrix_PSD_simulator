---
title: Tmatrix PSD Simulator
colorFrom: blue
colorTo: green
sdk: docker
app_port: 7860
pinned: false
---

# Local T-matrix PSD Explorer

Interactive polarimetric-radar PSD explorer backed by the legacy local
Fortran T-matrix code vendored in `engine/py_Tmatrix_Mueller`.

The app follows the same basic structure as
[myPSD](https://huggingface.co/spaces/snesbitt/myPSD): browser controls call a
JSON compute endpoint, the backend builds/caches T-matrix scatter tables, and
the frontend plots the PSD plus radar variables.

## Shareable Deployment

This app needs a Python backend and Fortran executables, so it cannot be shared
as a plain GitHub Pages static site. Use a Docker-capable host instead.

Recommended option: Hugging Face Spaces

1. Create a new Space.
2. Choose **Docker** as the SDK.
3. Connect or upload this GitHub repository.
4. The included `Dockerfile` installs `gfortran`, compiles the Fortran sources,
   and starts the app on port `7860`.

The public URL will look like:

```text
https://<username>-<space-name>.hf.space
```

Render, Railway, Fly.io, and other Docker hosts should also work. Set:

```text
PORT=7860
HOST=0.0.0.0
```

## Local Run

```bash
python3 -m tmatrix_web.backend.app.main
```

Open:

```text
http://127.0.0.1:7860
```

## Current Validation Status

Rain is the verified local path. It now mirrors the myPSD DSD handling:

- PSD curve and total number concentration use `D = 0.05-19.95 mm`.
- Rain radar-variable integration uses `D_max = 10 mm`.
- The legacy `Scatterer.radius` field is treated as equivalent diameter in mm,
  matching the bundled notebooks and the Fortran `DEQ` input.

Frozen rain, dry hail, spongy graupel, and spongy aggregate execute as legacy
diagnostic paths, but their dielectric/density assumptions still need
reproducible validation before scientific use.

See [tmatrix_web/README.md](tmatrix_web/README.md) for implementation details
and the comparison against myPSD/rustmatrix.
