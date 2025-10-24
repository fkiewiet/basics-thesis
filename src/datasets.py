# src/datasets.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Tuple, Iterable, Optional, Callable
import numpy as np

from .config import GridSpec
from .operators import Discretisation
from .loads import PointSource, RandomPointSource, build_load
from .solvers import GMRESOptions, gmres_solve

def _complex_to_channels(u: np.ndarray, shape: Tuple[int, ...]) -> np.ndarray:
    u = u.reshape(shape)
    return np.stack([u.real, u.imag], axis=0)  # (2,H,W)

def _scalar_field_channel(shape: Tuple[int, ...], value: float) -> np.ndarray:
    return np.full(shape, float(value), dtype=np.float32)

def build_direct_map(
    n: int,
    *,
    grid: GridSpec,
    disc: Discretisation,
    k_range: Tuple[float, float],
    gmres_opts: Optional[GMRESOptions] = None,
) -> Dict[str, np.ndarray]:
    """Return dict with (f,k)→u arrays:
       x = concat([Re f, Im f, k_field]); y = Re/Im u.
    """
    rng = np.random.default_rng(7)
    H, W = grid.shape
    xs, ys = [], []
    for i in range(n):
        k = float(rng.uniform(*k_range))
        # randomized point source (amp + phase)
        load = RandomPointSource(seed=rng.integers(0, 1_000_000))
        b = build_load(load, grid)                                      # flat RHS, complex128
        A = disc.assemble(grid, k)                                      # sparse matrix (Δ + k^2 I)
        res = gmres_solve(A, b, options=gmres_opts)                     # complex solution
        x = np.stack([
            _complex_to_channels(b, grid.shape)[0],
            _complex_to_channels(b, grid.shape)[1],
            _scalar_field_channel(grid.shape, k),
        ], axis=0)                                                      # (3,H,W)
        y = _complex_to_channels(res.solution, grid.shape)              # (2,H,W)
        xs.append(x.astype(np.float32)); ys.append(y.astype(np.float32))
    return {"x": np.stack(xs), "y": np.stack(ys)}                       # shapes: (N,3,H,W), (N,2,H,W)

def build_freq_transfer(
    n: int,
    *,
    grid: GridSpec,
    disc: Discretisation,
    k_low_range: Tuple[float, float],
    k_high_range: Tuple[float, float],
    gmres_opts: Optional[GMRESOptions] = None,
) -> Dict[str, np.ndarray]:
    """Return dict with (u_low,k_high)→u_high arrays:
       x = concat([Re u_low, Im u_low, k_high_field]); y = Re/Im u_high.
    """
    rng = np.random.default_rng(11)
    H, W = grid.shape
    xs, ys = [], []
    for i in range(n):
        k_low  = float(rng.uniform(*k_low_range))
        k_high = float(rng.uniform(*k_high_range))
        load = RandomPointSource(seed=rng.integers(0, 1_000_000))
        b = build_load(load, grid)                                      # same source for low/high
        A_low  = disc.assemble(grid, k_low)
        A_high = disc.assemble(grid, k_high)
        u_low  = gmres_solve(A_low,  b, options=gmres_opts).solution
        u_high = gmres_solve(A_high, b, options=gmres_opts).solution
        u_low_ch  = _complex_to_channels(u_low, grid.shape)
        u_high_ch = _complex_to_channels(u_high, grid.shape)
        x = np.stack([u_low_ch[0], u_low_ch[1], _scalar_field_channel(grid.shape, k_high)], axis=0)
        y = u_high_ch
        xs.append(x.astype(np.float32)); ys.append(y.astype(np.float32))
    return {"x": np.stack(xs), "y": np.stack(ys)}


SolverFn = Callable[..., object]  # expects (A, b, options?) -> object with .solution

def _complex_to_channels(u: np.ndarray, shape) -> np.ndarray:
    u = u.reshape(shape)
    return np.stack([u.real, u.imag], axis=0)

def _scalar_field_channel(shape, value: float) -> np.ndarray:
    return np.full(shape, float(value), dtype=np.float32)

def build_direct_map(
    n: int,
    *,
    grid: GridSpec,
    disc: Discretisation,
    k_range: Tuple[float, float],
    gmres_opts: Optional[GMRESOptions] = None,   # kept for API compat
    solver_fn: Optional[SolverFn] = None,        # <- NEW
    solver_options: Optional[object] = None,     # <- NEW (DirectOptions or GMRESOptions)
) -> Dict[str, np.ndarray]:
    """(f,k)→u pairs."""
    rng = np.random.default_rng(7)
    xs, ys = [], []
    for _ in range(n):
        k = float(rng.uniform(*k_range))
        load = RandomPointSource(seed=rng.integers(0, 1_000_000))
        b = build_load(load, grid)                # complex128 vector (flat)
        A = disc.assemble(grid, k)

        # choose backend
        if solver_fn is None:
            res = gmres_solve(A, b, options=gmres_opts)   # default = GMRES
        else:
            # pass solver_options if provided; otherwise default call
            res = solver_fn(A, b, solver_options) if solver_options is not None else solver_fn(A, b)

        u = res.solution
        x = np.stack([
            _complex_to_channels(b, grid.shape)[0],
            _complex_to_channels(b, grid.shape)[1],
            _scalar_field_channel(grid.shape, k),
        ], axis=0)
        y = _complex_to_channels(u, grid.shape)
        xs.append(x.astype(np.float32)); ys.append(y.astype(np.float32))
    return {"x": np.stack(xs), "y": np.stack(ys)}

def build_freq_transfer(
    n: int,
    *,
    grid: GridSpec,
    disc: Discretisation,
    k_low_range: Tuple[float, float],
    k_high_range: Tuple[float, float],
    gmres_opts: Optional[GMRESOptions] = None,   # kept for API compat
    solver_fn: Optional[SolverFn] = None,        # <- NEW
    solver_options: Optional[object] = None,     # <- NEW
) -> Dict[str, np.ndarray]:
    """(u_low, k_high)→u_high pairs."""
    rng = np.random.default_rng(11)
    xs, ys = [], []
    for _ in range(n):
        k_low  = float(rng.uniform(*k_low_range))
        k_high = float(rng.uniform(*k_high_range))
        load = RandomPointSource(seed=rng.integers(0, 1_000_000))
        b = build_load(load, grid)
        A_low  = disc.assemble(grid, k_low)
        A_high = disc.assemble(grid, k_high)

        if solver_fn is None:
            u_low  = gmres_solve(A_low,  b, options=gmres_opts).solution
            u_high = gmres_solve(A_high, b, options=gmres_opts).solution
        else:
            if solver_options is not None:
                u_low  = solver_fn(A_low,  b, solver_options).solution
                u_high = solver_fn(A_high, b, solver_options).solution
            else:
                u_low  = solver_fn(A_low,  b).solution
                u_high = solver_fn(A_high, b).solution

        u_low_ch  = _complex_to_channels(u_low, grid.shape)
        u_high_ch = _complex_to_channels(u_high, grid.shape)
        x = np.stack([u_low_ch[0], u_low_ch[1], _scalar_field_channel(grid.shape, k_high)], axis=0)
        y = u_high_ch
        xs.append(x.astype(np.float32)); ys.append(y.astype(np.float32))
    return {"x": np.stack(xs), "y": np.stack(ys)}