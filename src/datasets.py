# src/datasets.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple, Optional, Callable, Type

import numpy as np

from .config import GridSpec
from .operators import PMLConfig, assemble_operator, Discretisation  # Discretisation kept for legacy typing
from .loads import RandomPointSource, build_load
from .solvers import GMRESOptions, gmres_solve


# ------------------------- public typing -------------------------

# A solver that accepts (A, b, options?) and returns an object with `.solution` (GMRES-compatible).
SolverFn = Callable[..., object]


# ------------------------- small helpers -------------------------

def _complex_to_channels(u: np.ndarray, shape: Tuple[int, ...]) -> np.ndarray:
    """
    Reshape a flat complex array to (H,W) and split into (2,H,W) = [Re, Im].
    """
    u = np.asarray(u).reshape(shape)
    return np.stack([u.real, u.imag], axis=0)  # (2,H,W)


def _scalar_field_channel(shape: Tuple[int, ...], value: float) -> np.ndarray:
    """
    Produce a (H,W) scalar field channel filled with a constant (e.g., k).
    """
    return np.full(shape, float(value), dtype=np.float32)


# ------------------------- dataset builders -------------------------

def build_direct_map(
    n: int,
    *,
    grid: GridSpec,
    disc: Discretisation | None = None,  # kept for API-compat; not used by new assembly path
    k_range: Tuple[float, float],
    # solver selection
    gmres_opts: Optional[GMRESOptions] = None,       # used if solver_fn is None
    solver_fn: Optional[SolverFn] = None,            # plug in direct_solve(...) etc.
    solver_options: Optional[object] = None,         # DirectOptions or GMRESOptions for custom fn
    # physics / assembly
    pml: Optional[PMLConfig] = None,                 # <<< NEW: PML support
    dtype: np.dtype = np.complex128,
    # data variance
    seed: Optional[int] = 7,
) -> Dict[str, np.ndarray]:
    """
    Build (f, k) → u pairs for a Helmholtz problem on a rectangular grid.

    Returns
    -------
    dict with:
      - x: (N, 3, H, W)   = [Re f, Im f, k_field]
      - y: (N, 2, H, W)   = [Re u, Im u]
    """
    rng = np.random.default_rng(seed)
    H, W = grid.shape

    xs, ys = [], []

    for _ in range(n):
        # 1) sample frequency
        k = float(rng.uniform(*k_range))

        # 2) random point source (amplitude, phase, interior location)
        load = RandomPointSource(seed=int(rng.integers(0, 1_000_000)))
        b = build_load(load, grid)  # flat complex128 vector, length == grid.size

        # 3) assemble A(k) with optional PML
        A = assemble_operator(
            grid=grid,
            k=k,
            kind="helmholtz",
            dtype=dtype,
            pml=pml,  # <<< PML passed through
        )

        # 4) solve A u = b
        if solver_fn is None:
            res = gmres_solve(A, b, options=gmres_opts)  # default backend
        else:
            res = solver_fn(A, b, solver_options) if solver_options is not None else solver_fn(A, b)

        u = np.asarray(res.solution)

        # 5) pack channels
        f_ch = _complex_to_channels(b, grid.shape)          # (2,H,W)
        u_ch = _complex_to_channels(u, grid.shape)          # (2,H,W)
        x = np.stack([f_ch[0], f_ch[1], _scalar_field_channel(grid.shape, k)], axis=0)  # (3,H,W)

        xs.append(x.astype(np.float32, copy=False))
        ys.append(u_ch.astype(np.float32, copy=False))

    return {"x": np.stack(xs, axis=0), "y": np.stack(ys, axis=0)}  # shapes: (N,3,H,W), (N,2,H,W)


def build_freq_transfer(
    n: int,
    *,
    grid: GridSpec,
    disc: Discretisation | None = None,  # kept for API-compat; not used by new assembly path
    k_low_range: Tuple[float, float],
    k_high_range: Tuple[float, float],
    # solver selection
    gmres_opts: Optional[GMRESOptions] = None,       # used if solver_fn is None
    solver_fn: Optional[SolverFn] = None,            # plug in direct_solve(...) etc.
    solver_options: Optional[object] = None,         # DirectOptions or GMRESOptions for custom fn
    # physics / assembly
    pml: Optional[PMLConfig] = None,                 # <<< NEW: PML support
    dtype: np.dtype = np.complex128,
    # data variance
    seed: Optional[int] = 11,
) -> Dict[str, np.ndarray]:
    """
    Build (u_low, k_high) → u_high pairs using the SAME random source for low/high k.

    Returns
    -------
    dict with:
      - x: (N, 3, H, W)   = [Re u_low, Im u_low, k_high_field]
      - y: (N, 2, H, W)   = [Re u_high, Im u_high]
    """
    rng = np.random.default_rng(seed)
    H, W = grid.shape

    xs, ys = [], []

    for _ in range(n):
        # 1) sample low/high frequencies
        k_low = float(rng.uniform(*k_low_range))
        k_high = float(rng.uniform(*k_high_range))

        # 2) same randomized source for supervision
        load = RandomPointSource(seed=int(rng.integers(0, 1_000_000)))
        b = build_load(load, grid)

        # 3) assemble A_low, A_high with PML if requested
        A_low = assemble_operator(
            grid=grid,
            k=k_low,
            kind="helmholtz",
            dtype=dtype,
            pml=pml,
        )
        A_high = assemble_operator(
            grid=grid,
            k=k_high,
            kind="helmholtz",
            dtype=dtype,
            pml=pml,
        )

        # 4) solve
        if solver_fn is None:
            u_low = gmres_solve(A_low, b, options=gmres_opts).solution
            u_high = gmres_solve(A_high, b, options=gmres_opts).solution
        else:
            if solver_options is not None:
                u_low = solver_fn(A_low, b, solver_options).solution
                u_high = solver_fn(A_high, b, solver_options).solution
            else:
                u_low = solver_fn(A_low, b).solution
                u_high = solver_fn(A_high, b).solution

        # 5) pack channels
        u_low_ch = _complex_to_channels(u_low, grid.shape)      # (2,H,W)
        u_high_ch = _complex_to_channels(u_high, grid.shape)    # (2,H,W)
        x = np.stack([u_low_ch[0], u_low_ch[1], _scalar_field_channel(grid.shape, k_high)], axis=0)  # (3,H,W)

        xs.append(x.astype(np.float32, copy=False))
        ys.append(u_high_ch.astype(np.float32, copy=False))

    return {"x": np.stack(xs, axis=0), "y": np.stack(ys, axis=0)}
