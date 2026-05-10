"""Local T-matrix compute layer for the PSD explorer.

This module mirrors the myPSD backend split, but it uses the Fortran-backed
`py_Tmatrix_Mueller` package already present in this workspace instead of
`rustmatrix`.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import gamma
from pathlib import Path
import sys
import threading

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
ENGINE_ROOT = ROOT / "engine"
if ENGINE_ROOT.exists() and str(ENGINE_ROOT) not in sys.path:
    sys.path.insert(0, str(ENGINE_ROOT))

from py_Tmatrix_Mueller import radar, tmatrix_aux
from py_Tmatrix_Mueller.scatter import Scatterer


_TMATRIX_LOCK = threading.Lock()

_BANDS = {
    "S": tmatrix_aux.wl_S,
    "C": tmatrix_aux.wl_C,
    "X": tmatrix_aux.wl_X,
}

_D_PLOT = np.arange(0.05, 20.0, 0.05)


@dataclass(frozen=True)
class Hydrometeor:
    label: str
    dclass: int
    iopt: int
    d_max: float
    axis_ratio: str
    verified: bool = False
    executable: bool = True
    material_note: str = ""
    validation_note: str = ""


HYDROMETEORS = {
    "rain": Hydrometeor(
        label="Rain",
        dclass=1,
        iopt=1,
        d_max=10.0,
        axis_ratio="thurai_2007",
        verified=True,
        material_note="Liquid-water dielectric constant from EPSLON at the selected wavelength and temperature.",
        validation_note="Rain follows the myPSD assumptions most closely: oblate drops, liquid water at 10 C, D_max=10 mm, and the Testud normalized-gamma PSD. Remaining differences can come from rustmatrix/pytmatrix versus this legacy Fortran Mueller implementation.",
    ),
    "frozen_rain": Hydrometeor(
        label="Frozen rain",
        dclass=2,
        iopt=2,
        d_max=10.0,
        axis_ratio="sphere",
        material_note="Legacy executable still uses the liquid-water EPSLON path unless py_tmat.f90 is recompiled with the non-rain dielectric overrides.",
        validation_note="Mueller file handling is fixed, but the physical dielectric assumption is not yet validated.",
    ),
    "dry_hail": Hydrometeor(
        label="Dry hail",
        dclass=3,
        iopt=3,
        d_max=40.0,
        axis_ratio="near_sphere",
        material_note="Source comments indicate dry-ice epsilon should be about 3.168351 + 0.02492i; current bundled executable cannot be confirmed to use it.",
        validation_note="Executable diagnostic only until the Fortran dielectric branch is rebuilt and compared against a reference case.",
    ),
    "spongy_graupel": Hydrometeor(
        label="Spongy graupel",
        dclass=7,
        iopt=2,
        d_max=15.0,
        axis_ratio="near_sphere",
        material_note="Original wrapper maps this through the frozen-raindrop/spongy-graupel T-matrix option; effective-medium details are not implemented in Python.",
        validation_note="Executable diagnostic only; density and ice-air/water fraction must be specified before scientific use.",
    ),
    "snow_aggregate": Hydrometeor(
        label="Spongy aggregate",
        dclass=11,
        iopt=8,
        d_max=15.0,
        axis_ratio="oblate_aggregate",
        material_note="Original comments list fixed aggregate epsilon examples by density, but the bundled executable is not confirmed to apply them.",
        validation_note="Executable diagnostic only; unlike myPSD, this does not compute per-diameter Maxwell-Garnett snow refractive index.",
    ),
}


def _axis_ratio(kind: str, diameter_mm: np.ndarray) -> np.ndarray:
    if kind == "thurai_2007":
        return tmatrix_aux.dsr_thurai_2007(diameter_mm)
    if kind == "near_sphere":
        return np.full_like(diameter_mm, 0.99, dtype=float)
    if kind == "oblate_aggregate":
        return np.full_like(diameter_mm, 1.25, dtype=float)
    return np.ones_like(diameter_mm, dtype=float)


@lru_cache(maxsize=32)
def _scatter_table(
    band: str,
    hydro_key: str,
    d_min: float,
    d_max: float,
    d_step: float,
    canting_std_deg: float,
) -> Scatterer:
    hydro = HYDROMETEORS[hydro_key]
    diameter = np.arange(d_min, d_max + 0.5 * d_step, d_step)
    scatterer = Scatterer(
        # The legacy wrapper calls this field `radius`, but its examples pass
        # equivalent diameter in mm and write it to Fortran as DEQ in cm.
        radius=diameter,
        axis_ratio=_axis_ratio(hydro.axis_ratio, diameter),
        wavelength=_BANDS[band],
        tempt=10.0,
        IB=8,
        IOPT=hydro.iopt,
        NM=4,
        NRANK=8,
        anginc=90.0,
        dclass=hydro.dclass,
        distyp=2,
        dmin=float(diameter[0]),
        dmax=float(diameter[-1]),
        dpar1=0.0,
        dpar2=max(float(canting_std_deg), 0.005),
        dpar3=0.0,
    )

    # The legacy Fortran wrapper writes fixed input/output filenames inside
    # py_Tmatrix_Mueller/fortran_tm, so table builds must be serialized.
    with _TMATRIX_LOCK:
        scatterer.get_table()
    return scatterer


def normalized_gamma(diameter_mm: np.ndarray, dm: float, log_nw: float, mu: float) -> np.ndarray:
    f_mu = (6.0 / 4.0**4) * ((4.0 + mu) ** (mu + 4.0)) / gamma(mu + 4.0)
    nd = (
        10.0**log_nw
        * f_mu
        * (diameter_mm / dm) ** mu
        * np.exp(-(4.0 + mu) * (diameter_mm / dm))
    )
    return np.nan_to_num(nd, nan=0.0, posinf=0.0, neginf=0.0)


def compute(
    dm: float,
    log_nw: float,
    mu: float,
    band: str,
    hydro_key: str,
    canting_std_deg: float,
    d_step: float = 0.1,
) -> dict:
    if band not in _BANDS:
        raise ValueError(f"Unsupported band: {band}")
    if hydro_key not in HYDROMETEORS:
        raise ValueError(f"Unsupported hydrometeor: {hydro_key}")

    hydro = HYDROMETEORS[hydro_key]
    d_min = 0.1
    d_max = hydro.d_max
    diameter = np.arange(d_min, d_max + 0.5 * d_step, d_step)
    nd = normalized_gamma(diameter, dm, log_nw, mu)
    nd_plot = normalized_gamma(_D_PLOT, dm, log_nw, mu)
    scatterer = _scatter_table(
        band,
        hydro_key,
        d_min,
        d_max,
        d_step,
        round(float(canting_std_deg), 3),
    )

    zh = float(radar.calc_ZH(scatterer, nd, diameter, d_step))
    zv = float(radar.calc_ZV(scatterer, nd, diameter, d_step))
    zdr = float(radar.calc_ZDR(scatterer, nd, diameter, d_step))
    kdp = float(radar.calc_KDP(scatterer, nd, diameter, d_step))
    rhohv = float(radar.calc_RHOHV(scatterer, nd, diameter, d_step))
    ah = float(radar.calc_ATTH(scatterer, nd, diameter, d_step))
    adr = float(radar.calc_DATT(scatterer, nd, diameter, d_step))
    nt = float(np.trapz(nd_plot, _D_PLOT))

    assumption_bullets = [
        "T-matrix table is produced by the local Fortran executables in py_Tmatrix_Mueller/fortran_tm.",
        f"Legacy species switches: dclass={hydro.dclass}, IOPT={hydro.iopt}.",
        hydro.material_note,
        hydro.validation_note,
        f"Axis-ratio model: {hydro.axis_ratio}; PSD integration D = {d_min:.1f}-{d_max:.1f} mm.",
        "PSD form follows the Testud normalized-gamma convention used by myPSD.",
    ]
    if hydro_key == "rain":
        assumption_bullets.append(
            "Rain reference: Beard and Chuang (1987) equilibrium drop shape; Testud et al. (2001) normalized-gamma PSD; myPSD uses rustmatrix.refractive.m_w_10C."
        )

    return {
        "metrics": {
            "zh_dbz": zh,
            "zv_dbz": zv,
            "zdr_db": zdr,
            "kdp_deg_per_km": kdp,
            "rho_hv": rhohv,
            "ah_db_per_km": ah,
            "adr_db_per_km": adr,
            "nt_per_m3": nt,
        },
        "nd": {
            "d_mm": _D_PLOT.tolist(),
            "n_d": nd_plot.tolist(),
        },
        "assumptions": {
            "title": hydro.label,
            "verified": hydro.verified,
            "executable": hydro.executable,
            "bullets": assumption_bullets,
        },
    }
