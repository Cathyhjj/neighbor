"""Batch neighbor searches with explicit periodic boundary handling.

pymatgen is optional. Without it, use SciPy for finite clusters and ASE for
periodic structures. All backends return directed edges and actual image
distances, including nonzero periodic images of an atom itself.
"""

import numpy as np
from ase.cell import Cell
from ase.neighborlist import neighbor_list
from scipy.spatial import cKDTree


BACKENDS = {"auto", "pymatgen", "scipy", "ase"}


def find_neighbor_edges(atoms, radius, pbc, backend="auto"):
    """Return ``(i, j, images, distances, backend_used)`` within ``radius``.

    Coordinates are Cartesian in Angstrom. The search includes the radius
    boundary; callers apply their own shell intervals. Only zero-image self
    edges are removed. Coincident distinct sites are rejected.
    """
    if backend not in BACKENDS:
        raise ValueError("neighbor_backend must be auto, pymatgen, scipy or ase")
    radius = float(radius)
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("Neighbor search radius must be finite and positive")
    xyz = np.ascontiguousarray(atoms.positions, dtype=float)
    pbc = np.asarray(pbc, dtype=bool)
    if not np.isfinite(xyz).all():
        raise ValueError("Atomic positions must be finite")
    if pbc.any():
        cell = np.asarray(atoms.cell, dtype=float)
        if (not np.isfinite(cell).all()
                or np.linalg.matrix_rank(cell[pbc]) != np.count_nonzero(pbc)):
            raise ValueError("Periodic directions require independent, nonzero cell vectors")
        # Nonperiodic vectors do not define image translations. They may be
        # dependent on the periodic vectors (e.g. a chain with a redundant box
        # axis), so discard them before completing such a lattice.
        if np.linalg.matrix_rank(cell) < np.count_nonzero(np.any(cell != 0, axis=1)):
            cell = cell.copy()
            cell[~pbc] = 0
        lattice = np.ascontiguousarray(Cell(cell).complete(), dtype=float)
        if not np.isfinite(lattice).all() or np.linalg.matrix_rank(lattice) != 3:
            raise ValueError("Could not construct a nonsingular cell for periodic neighbor search")
    else:
        lattice = np.eye(3, dtype=float)

    search = None
    if backend in {"auto", "pymatgen"}:
        try:
            from pymatgen.optimization.neighbors import find_points_in_spheres
        except ImportError as error:
            if backend == "pymatgen":
                raise ImportError("Install pymatgen to use neighbor_backend='pymatgen'") from error
        else:
            search = find_points_in_spheres
    used = "pymatgen" if search is not None else backend
    if used == "auto":
        used = "ase" if pbc.any() else "scipy"
    if used == "scipy" and pbc.any():
        raise ValueError("The scipy backend supports finite clusters only; use auto, ase or pymatgen for PBC")

    if len(atoms) == 0:
        return (np.empty(0, dtype=int), np.empty(0, dtype=int),
                np.empty((0, 3), dtype=int), np.empty(0), used)

    # ASE uses a strict cutoff, so include the requested endpoint in candidates.
    search_radius = radius + max(1e-10, radius * 1e-12)
    if used == "pymatgen":
        i, j, images, distances = search(
            xyz, xyz, search_radius,
            pbc=np.ascontiguousarray(pbc, dtype=np.int64), lattice=lattice,
        )
    elif used == "scipy":
        pairs = cKDTree(xyz).query_pairs(search_radius, output_type="ndarray")
        first, second = pairs.T
        lengths = np.linalg.norm(xyz[second] - xyz[first], axis=1)
        i, j = np.r_[first, second], np.r_[second, first]
        images = np.zeros((len(i), 3), dtype=int)
        distances = np.r_[lengths, lengths]
    else:
        work = atoms.copy()
        work.set_pbc(pbc)
        if not pbc.any():
            # Ordinary XYZ files have no cell. Bound the points so ASE's
            # spatial bins do not degenerate into an all-pairs calculation.
            work.positions -= xyz.min(axis=0) - search_radius
            work.set_cell(np.ptp(xyz, axis=0) + 2 * search_radius)
        else:
            work.set_cell(lattice)
        i, j, images, distances = neighbor_list("ijSd", work, search_radius)

    i, j = np.asarray(i, dtype=int), np.asarray(j, dtype=int)
    images = np.asarray(images, dtype=int).reshape(-1, 3)
    distances = np.asarray(distances, dtype=float)
    self_edge = (i == j) & np.all(images == 0, axis=1)
    keep = ~self_edge & (distances <= radius)
    if np.any(distances[keep] <= 1e-12):
        raise ValueError("Coincident atomic sites found; CN is undefined for overlapping sites")
    return i[keep], j[keep], images[keep], distances[keep], used
