Neighbor: This calculates the CNs information for a cluster
-----------------------------------------------------------------------------

**Try out in online Jupyter Notebook**: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Cathyhjj/neighbor.git/main?labpath=20240806_simple_example.ipynb) or [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1dvgHGLjOLWhGTAJJKKPIf3th4PG9YETN?usp=sharing)

author:    Juanjuan Huang & Shelly D. Kelly

email:     juanjuan.huang@anl.gov

copyright: Copyright © 2024, UChicago Argonne, LLC

## Coordination numbers

Existing notebook calls to `ClusterNeighbor()`, `get_CN()`, and
`get_CN_summary_all()` still work. The CN routines share a batch neighbor table
and one shell assignment, so shell means, per-atom counts, and distributions
refer to the same neighbors. Counts are directed: `Ce-O` means O neighbors per
Ce atom. All central atoms, including those with zero neighbors, enter the mean.

Use the dependencies in `environment.yml`. For the optional pymatgen search
backend, install it in the same Python environment:

```sh
python -m pip install pymatgen
```

```python
from neighbor_core import ClusterNeighbor

cluster = ClusterNeighbor(neighbor_backend="auto", boundary_mode="cluster")
cluster.load_xyz(path="test_Pt/testPt.xyz")
cluster.get_CN_summary_all(tolerance=0.1, bond_range=5, printit=False)
print(cluster.CN["Pt-Pt"])  # Mean distance in Å -> mean CN
print(cluster.neighbor_backend_used)
```

`neighbor_backend="auto"` prefers pymatgen's `find_points_in_spheres`. When
pymatgen is unavailable, it uses SciPy's `cKDTree` for finite clusters and ASE
for periodic structures. You can explicitly select `"pymatgen"`, `"scipy"`, or
`"ase"`; SciPy supports finite clusters only. Search backends change how neighbors
are found, not the definition of radial-shell CN.

Choose the boundary conditions for the structure being analyzed:

| `boundary_mode` | Behavior |
| --- | --- |
| `"auto"` (default) | Respect `atoms.pbc`, including partial periodicity. |
| `"cluster"` | Disable periodicity for a finite cluster. |
| `"periodic"` | Enable periodicity in all three directions; requires a valid cell. |

Periodic calculations retain image offsets and count nonzero images of the
central atom when appropriate. `expand_to_sphere()` and `shrink_cluster_size()`
return finite clusters with PBC disabled; use `boundary_mode="auto"` or
`"cluster"` for these cut structures.

## Defining radial shells

The default groups sorted distances into shells wherever the adjacent distance
gap exceeds `tolerance` (Å). This preserves the existing grouping rule, including
its tendency to join a chain of nearby distances. For distorted clusters or a
series of particle sizes, explicit `shell_edges` provide consistent intervals:

```python
edges = {"Pt-Pt": [0, 3.3, 4.2, 5.0]}  # Illustrative bounds in Å
cluster.get_CN_summary_all(bond_range=5, shell_edges=edges, printit=False)

first_shell = cluster.CN_shells["Pt-Pt"][0]
print(first_shell["shell_id"], first_shell["average_CN"])
print(first_shell["atom_indices"], first_shell["atom_CN"])
print(first_shell["distribution"])
print(first_shell["distance"], first_shell["distance_std"])
```

Each interval is **[lower, upper)**: the lower boundary is included and the
upper boundary is excluded. Edges must be finite, nonnegative, strictly
increasing, and end at or below `bond_range`. All radial-shell calculations
exclude distances at or above `bond_range`. Choose boundaries appropriate to
your structure; the example is not a universal bond criterion.

A sequence applies to every requested pair. A mapping accepts string keys such
as `"Ce-O"` or tuple keys such as `("Ce", "O")`; the reverse pair reuses the same
edges unless it has its own entry. Pairs missing from the mapping use the default
gap grouping.

`CN_shells[pair]` contains one record per shell, with 1-based `shell_id`, bounds,
method, mean/minimum/maximum distance, population `distance_std`, `average_CN`,
and NumPy arrays `atom_indices`, `atom_CN`, and `distances`. `distribution` maps
CN to the number of central atoms. Explicit bins retain their IDs even when
empty: their distance statistics are `None`, and their CN is zero when central
atoms exist. The legacy distance-keyed `CN` and the shell summaries in
`CN_summary` contain occupied shells only. If the central element is absent,
its mean CN is undefined (`None`)
and its arrays/distribution are empty; calculations over all pairs of an empty
structure return empty results. Single-pair queries on an empty structure need
explicit element names.

To obtain a distribution for an existing shell, pass the same calculation
parameters, including `shell_edges`:

```python
distribution = cluster.get_CN_around_distance(
    "Pt", "Pt", shell=1, bond_range=5, shell_edges=edges, printit=False
)
```

`shell=1` uses the exact first-shell members. Alternatively,
`target_distance=2.77, tolerance=0.05` explicitly requests the inclusive distance
window `[2.72, 2.82]` Å, still restricted to distances below `bond_range`. A target
window is a separate selection and can differ from a grouped shell; do not pass
both `shell` and `target_distance`.

`get_CN_for_atom_index()` uses an inclusive pair-distance `cutoff`; when omitted,
it uses the sum of the two atoms' covalent radii. Its `tolerance` argument remains
accepted for compatibility and does not broaden the cutoff.
`remove_under_coordinated_atoms()` uses the same cutoff convention and now
includes isolated (CN=0) atoms by default. Pass `include_isolated=False` to keep
them. Removal is one pass, rather than repeated pruning, and `self_apply=False`
returns the reduced structure without replacing the loaded one.

Changes to positions, species, cell, PBC, masses, or search options invalidate
cached results on the next CN calculation, including direct edits to
`cluster.atoms`. Call `refresh_atoms()` to clear derived results immediately.

`write_to_excel()` preserves the original columns and adds shell ID, method,
bounds, and distance standard deviation. It exports the current summaries
(or calculated shells if no summaries were requested), including empty explicit
bins with zero CN and blank distance statistics. For `method="gap"`, bounds
are the minimum and maximum observed distances, not fixed shell cutoffs.
For `method="edges"`, they are the supplied half-open interval boundaries.

## Validation

From the repository root, with the project dependencies installed:

```sh
python -m unittest discover -s tests -v
```

The tests cover finite and periodic coordination, shell membership, empty bins,
cache invalidation, and backend agreement. Pymatgen-specific tests are skipped
when that optional dependency is unavailable.
