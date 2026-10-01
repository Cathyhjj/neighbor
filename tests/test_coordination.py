"""Numerical regressions for radial coordination numbers.

Run with ``python -m unittest discover -s tests -v`` from the repository root.
Tests use exact geometric counts and an independent finite-cluster reference.
"""

import builtins
import contextlib
import io
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io import read

import neighbor_core
from neighbor_core import ClusterNeighbor


try:
    from pymatgen.optimization.neighbors import find_points_in_spheres
except ImportError:
    HAVE_PYMATGEN = False
else:
    HAVE_PYMATGEN = True


ROOT = Path(__file__).resolve().parents[1]


def make_cluster(atoms, backend="auto", boundary_mode="auto"):
    cluster = ClusterNeighbor(neighbor_backend=backend, boundary_mode=boundary_mode)
    cluster.load_xyz(from_file=False, atom_object=atoms)
    return cluster


def quiet(function, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return function(*args, **kwargs)


def finite_shell_counts(atoms, center_symbol, neighbor_symbol, edges):
    """Small brute-force oracle: every center, including zero-degree atoms."""
    symbols = atoms.get_chemical_symbols()
    centers = [i for i, symbol in enumerate(symbols) if symbol == center_symbol]
    result = np.zeros((len(edges) - 1, len(centers)), dtype=int)
    for column, i in enumerate(centers):
        for j, symbol in enumerate(symbols):
            if i == j or symbol != neighbor_symbol:
                continue
            distance = np.linalg.norm(atoms.positions[j] - atoms.positions[i])
            for shell, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
                if lo <= distance < hi:
                    result[shell, column] += 1
    return centers, result


class PeriodicCoordinationTests(unittest.TestCase):
    def test_primitive_fcc_retains_periodic_self_images(self):
        atoms = bulk("Pt", "fcc", a=3.92)
        self.assertEqual(len(atoms), 1)
        for backend in ("auto", "ase"):
            with self.subTest(backend=backend):
                cluster = make_cluster(atoms.copy(), backend)
                values = cluster.get_CN("Pt", "Pt", bond_range=5, printit=False)["Pt-Pt"]
                np.testing.assert_allclose(list(values), [3.92 / np.sqrt(2), 3.92, 3.92 * np.sqrt(1.5)])
                self.assertEqual(list(values.values()), [12, 6, 24])
                self.assertEqual(cluster.get_CN_for_atom_index(0, "Pt-Pt", cutoff=3.2), 12)
                self.assertEqual(cluster.get_CN_around_distance("Pt", "Pt", shell=1,
                                                               bond_range=5, printit=False), {12: 1})

    def test_repository_pt_cif_first_shell_is_twelve(self):
        atoms = read(ROOT / "structure_examples" / "Ptfoil_243678.cif")
        self.assertTrue(np.all(atoms.pbc))
        cluster = make_cluster(atoms)
        result = cluster.get_CN("Pt", "Pt", bond_range=3.2, printit=False)
        self.assertEqual(list(result["Pt-Pt"].values()), [12])
        self.assertEqual(cluster.get_CN_around_distance("Pt", "Pt", shell=1,
                                                       bond_range=3.2, printit=False), {12: len(atoms)})

    def test_boundary_mode_overrides_atoms_pbc(self):
        primitive = bulk("Pt", "fcc", a=3.92)
        cluster = make_cluster(primitive.copy(), boundary_mode="cluster")
        self.assertEqual(cluster.get_CN("Pt", "Pt", bond_range=3.2, printit=False)["Pt-Pt"], {})
        primitive.pbc = False
        periodic = make_cluster(primitive, boundary_mode="periodic")
        self.assertEqual(list(periodic.get_CN("Pt", "Pt", bond_range=3.2,
                                             printit=False)["Pt-Pt"].values()), [12])

    def test_partial_pbc_has_only_chain_neighbors(self):
        atoms = Atoms("Cu", positions=[[0, 0, 0]], cell=[2, 2, 2], pbc=[True, False, False])
        for backend in ("auto", "ase"):
            with self.subTest(backend=backend):
                cluster = make_cluster(atoms.copy(), backend)
                self.assertEqual(cluster.get_CN("Cu", "Cu", bond_range=2.5,
                                                printit=False)["Cu-Cu"], {2.0: 2.0})

    def test_partial_pbc_ignores_dependent_nonperiodic_cell_vectors(self):
        cases = [([[2, 0, 0], [4, 0, 0], [0, 0, 0]], [True, False, False], 2),
                 ([[2, 0, 0], [0, 2, 0], [2, 2, 0]], [True, True, False], 4)]
        for cell, pbc, expected in cases:
            for backend in ("auto", "ase"):
                with self.subTest(pbc=pbc, backend=backend):
                    cluster = make_cluster(Atoms("Cu", positions=[[0, 0, 0]], cell=cell, pbc=pbc), backend)
                    self.assertEqual(cluster.get_CN_for_atom_index(0, "Cu-Cu", cutoff=2.1), expected)

    def test_cut_cluster_outputs_disable_pbc(self):
        cluster = make_cluster(bulk("Pt", "fcc", a=3.92, cubic=True))
        sphere = cluster.expand_to_sphere(target_diameter=10)
        self.assertGreater(len(sphere), 0)
        self.assertFalse(np.any(sphere.pbc))
        self.assertTrue(np.all(cluster.atoms.pbc))
        cut = cluster.shrink_cluster_size(new_radius=2, center_atom_index=0)
        self.assertFalse(np.any(cut.pbc))
        cluster.expand_to_sphere(target_diameter=10, self_apply=True)
        self.assertFalse(np.any(cluster.atoms.pbc))


class ShellCoordinationTests(unittest.TestCase):
    def test_gap_group_membership_is_shared_by_all_statistics(self):
        atoms = Atoms("CuO5", positions=[[0, 0, 0], [2, 0, 0], [-2.009, 0, 0],
                                         [0, 2.018, 0], [0, -2.027, 0], [0, 0, 2.036]])
        cluster = make_cluster(atoms)
        result = cluster.get_CN("Cu", "O", tolerance=.01, bond_range=3, printit=False)
        self.assertEqual(list(result["Cu-O"].values()), [5])
        self.assertEqual(cluster.get_CN_around_distance("Cu", "O", shell=1,
                                                       tolerance=.01, bond_range=3, printit=False), {5: 1})
        summary = cluster.get_CN_summary_all(tolerance=.01, bond_range=3, printit=False)
        self.assertEqual(next(iter(summary["Cu-O"].values()))["average_CN"], 5)
        np.testing.assert_array_equal(cluster.CN_shells["Cu-O"][0]["atom_CN"], [5])

    def test_average_includes_zero_degree_centers_and_directed_denominator(self):
        atoms = Atoms("Cu2O", positions=[[0, 0, 0], [20, 0, 0], [2, 0, 0]])
        cluster = make_cluster(atoms)
        result = cluster.get_CN_all(bond_range=3, printit=False)
        self.assertIsInstance(result, dict)
        self.assertEqual(result["Cu-O"], {2.0: .5})
        self.assertEqual(result["O-Cu"], {2.0: 1.0})
        self.assertEqual(result["Cu-Cu"], {})
        self.assertEqual(result["O-O"], {})
        self.assertEqual(cluster.get_CN_around_distance("Cu", "O", shell=1,
                                                       bond_range=3, printit=False), {0: 1, 1: 1})

    def test_explicit_edges_preserve_empty_shell_ids_and_half_open_boundaries(self):
        atoms = Atoms("CuO3", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0], [4, 0, 0]])
        cluster = make_cluster(atoms)
        edges = [0, 1, 2, 3, 4]
        result = cluster.get_CN("Cu", "O", bond_range=4, printit=False, shell_edges=edges)
        shells = cluster.CN_shells["Cu-O"]
        self.assertEqual([shell["shell_id"] for shell in shells], [1, 2, 3, 4])
        self.assertEqual([shell["distance"] for shell in shells], [None, 1.0, 2.0, None])
        self.assertEqual([shell["average_CN"] for shell in shells], [0, 1, 1, 0])
        self.assertEqual(result["Cu-O"], {1.0: 1.0, 2.0: 1.0})
        for shell, lo, hi in zip(shells, edges[:-1], edges[1:]):
            self.assertEqual(shell["lower_bound"], lo)
            self.assertEqual(shell["upper_bound"], hi)
            self.assertEqual(list(shell["atom_indices"]), [0])
            self.assertEqual(shell["distribution"], {int(shell["average_CN"]): 1})
        self.assertEqual(cluster.get_CN_around_distance("Cu", "O", shell=4,
                                                       bond_range=4, printit=False,
                                                       shell_edges=edges), {0: 1})

    def test_pair_specific_edges_and_distance_spread(self):
        atoms = Atoms("Cu2O3", positions=[[0, 0, 0], [10, 0, 0], [1, 0, 0],
                                          [-2, 0, 0], [12, 0, 0]])
        cluster = make_cluster(atoms)
        edges = {"Cu-O": [0, 3, 4], "O-Cu": [0, 1.5, 3, 4]}
        cluster.get_CN_all(bond_range=4, shell_edges=edges, printit=False)
        shell = cluster.CN_shells["Cu-O"][0]
        self.assertAlmostEqual(shell["distance"], 5 / 3)
        self.assertAlmostEqual(shell["distance_std"], np.std([1, 2, 2]))
        self.assertEqual(shell["average_CN"], 1.5)
        np.testing.assert_array_equal(shell["atom_CN"], [2, 1])
        self.assertEqual([s["average_CN"] for s in cluster.CN_shells["O-Cu"]], [1 / 3, 2 / 3, 0])

    def test_reverse_pair_reuses_shell_edges(self):
        cluster = make_cluster(Atoms("CuO2", positions=[[0, 0, 0], [1.5, 0, 0], [-1.5, 0, 0]]))
        cluster.get_CN_all(bond_range=3, shell_edges={"Cu-O": [0, 1, 2, 3]}, printit=False)
        self.assertEqual([s["average_CN"] for s in cluster.CN_shells["Cu-O"]], [0, 2, 0])
        self.assertEqual([s["average_CN"] for s in cluster.CN_shells["O-Cu"]], [0, 1, 0])
        self.assertEqual([s["shell_id"] for s in cluster.CN_shells["O-Cu"]], [1, 2, 3])

    def test_mutating_shell_edges_in_place_reassigns_membership(self):
        cluster = make_cluster(Atoms("CuO3", positions=[[0, 0, 0], [1, 0, 0],
                                                      [2, 0, 0], [3, 0, 0]]))
        edges = np.array([0, 1.5, 4])
        cluster.get_CN("Cu", "O", bond_range=4, shell_edges=edges, printit=False)
        self.assertEqual([s["average_CN"] for s in cluster.CN_shells["Cu-O"]], [1, 2])
        edges[1] = 2.5
        cluster.get_CN("Cu", "O", bond_range=4, shell_edges=edges, printit=False)
        self.assertEqual([s["average_CN"] for s in cluster.CN_shells["Cu-O"]], [2, 1])

    def test_explicit_distance_window_stays_distinct_from_gap_shell(self):
        cluster = make_cluster(Atoms("CuO5", positions=[[0, 0, 0], [2, 0, 0],
                               [-2.009, 0, 0], [0, 2.018, 0], [0, -2.027, 0], [0, 0, 2.036]]))
        distribution = cluster.get_CN_around_distance("Cu", "O", target_distance=2.018,
                                                       tolerance=.01, bond_range=3, printit=False)
        self.assertEqual(distribution, {3: 1})

    def test_explicit_distance_window_is_clipped_by_search_radius(self):
        cluster = make_cluster(Atoms("CuO3", positions=[[0, 0, 0], [1.75, 0, 0],
                                                      [2, 0, 0], [2.25, 0, 0]]))
        distribution = cluster.get_CN_around_distance("Cu", "O", target_distance=2,
                                                       tolerance=.3, bond_range=2, printit=False)
        self.assertEqual(distribution, {1: 1})

    def test_bad_shell_edges_are_rejected(self):
        cluster = make_cluster(Atoms("Cu2", positions=[[0, 0, 0], [2, 0, 0]]))
        for edges in ([], [0], [0, 2, 2], [2, 1], [-1, 2], [0, float("nan")], [0, 6]):
            with self.subTest(edges=edges), self.assertRaises((TypeError, ValueError)):
                cluster.get_CN("Cu", "Cu", bond_range=5, printit=False, shell_edges=edges)


class CacheAndEdgeCaseTests(unittest.TestCase):
    def test_all_pair_summaries_share_one_physical_neighbor_search(self):
        cluster = make_cluster(Atoms("Cu2O2", positions=[[0, 0, 0], [2, 0, 0],
                                                       [0, 2, 0], [3, 3, 0]]))
        with patch("neighbor_core.find_neighbor_edges", wraps=neighbor_core.find_neighbor_edges) as search:
            cluster.get_CN_summary_all(bond_range=5, printit=False)
            self.assertEqual(search.call_count, 1)
            cluster.get_CN_summary_all(bond_range=3, tolerance=.1, printit=False)
            self.assertEqual(search.call_count, 1)
            cluster.atoms.positions[0, 0] += .1
            cluster.get_CN_summary_all(bond_range=3, tolerance=.1, printit=False)
            self.assertEqual(search.call_count, 2)

    def test_reducing_cutoff_clears_previous_pair_and_summary_results(self):
        cluster = make_cluster(Atoms("Cu3", positions=[[0, 0, 0], [2, 0, 0], [4, 0, 0]]))
        cluster.get_CN_summary_all(bond_range=5, printit=False)
        summary = cluster.get_CN_summary_all(bond_range=3, printit=False)
        self.assertEqual(set(cluster.CN["Cu-Cu"]), {2})
        self.assertEqual(set(summary["Cu-Cu"]), {2})
        cluster.get_CN("Cu", "Cu", bond_range=1, printit=False)
        self.assertEqual(cluster.CN["Cu-Cu"], {})
        self.assertEqual(cluster.CN_shells["Cu-Cu"], [])

    def test_direct_positions_and_species_changes_invalidate_results(self):
        cluster = make_cluster(Atoms("Cu2", positions=[[0, 0, 0], [2, 0, 0]]))
        self.assertEqual(cluster.get_CN("Cu", "Cu", bond_range=3, printit=False)["Cu-Cu"], {2: 1})
        cluster.atoms.positions[1] = [4, 0, 0]
        self.assertEqual(cluster.get_CN("Cu", "Cu", bond_range=3, printit=False)["Cu-Cu"], {})
        cluster.atoms.positions[1] = [2, 0, 0]
        cluster.atoms[1].symbol = "O"
        result = cluster.get_CN_all(bond_range=3, printit=False)
        self.assertEqual(result["Cu-Cu"], {})
        self.assertEqual(result["Cu-O"], {2: 1})
        self.assertEqual(cluster.elements_num, {"Cu": 1, "O": 1})

    def test_direct_pbc_and_cell_changes_invalidate_results(self):
        atoms = Atoms("Cu", positions=[[0, 0, 0]], cell=[2, 2, 2], pbc=True)
        cluster = make_cluster(atoms)
        self.assertEqual(cluster.get_CN_for_atom_index(0, "Cu-Cu", cutoff=2.5), 6)
        cluster.atoms.cell = [3, 3, 3]
        self.assertEqual(cluster.get_CN_for_atom_index(0, "Cu-Cu", cutoff=2.5), 0)
        cluster.atoms.cell = [2, 2, 2]
        cluster.atoms.pbc = False
        self.assertEqual(cluster.get_CN_for_atom_index(0, "Cu-Cu", cutoff=2.5), 0)

    def test_backend_and_boundary_mode_edits_invalidate_search(self):
        finite = make_cluster(Atoms("Cu2", positions=[[0, 0, 0], [2, 0, 0]]), "scipy")
        finite.get_CN_all(bond_range=3, printit=False)
        self.assertEqual(finite.neighbor_backend_used, "scipy")
        finite.neighbor_backend = "ase"
        self.assertEqual(finite.get_CN_all(bond_range=3, printit=False)["Cu-Cu"], {2: 1})
        self.assertEqual(finite.neighbor_backend_used, "ase")

        periodic = make_cluster(Atoms("Cu", positions=[[0, 0, 0]], cell=[2, 2, 2],
                                      pbc=[True, False, False]))
        self.assertEqual(periodic.get_CN_for_atom_index(0, "Cu-Cu", cutoff=2.5), 2)
        periodic.boundary_mode = "cluster"
        self.assertEqual(periodic.get_CN_for_atom_index(0, "Cu-Cu", cutoff=2.5), 0)
        periodic.boundary_mode = "periodic"
        self.assertEqual(periodic.get_CN_for_atom_index(0, "Cu-Cu", cutoff=2.5), 6)

    def test_single_atom_cutoff_is_inclusive_bulk_search_is_exclusive(self):
        cluster = make_cluster(Atoms("Cu2", positions=[[0, 0, 0], [2, 0, 0]]))
        self.assertEqual(cluster.get_CN("Cu", "Cu", bond_range=2, printit=False)["Cu-Cu"], {})
        self.assertEqual(cluster.get_CN_for_atom_index(0, "Cu-Cu", cutoff=2), 1)
        self.assertEqual(cluster.get_CN_for_atom_index(1, "Cu-Cu", cutoff=2), 1)

    def test_invalid_positions_cell_and_overlapping_sites_reject(self):
        backends = ["scipy", "ase"] + (["pymatgen"] if HAVE_PYMATGEN else [])
        for backend in backends:
            for positions in ([[0, 0, 0], [0, 0, 0]], [[0, 0, 0], [float("nan"), 0, 0]]):
                with self.subTest(backend=backend, positions=positions):
                    cluster = make_cluster(Atoms("Cu2", positions=positions), backend)
                    with self.assertRaises(ValueError):
                        cluster.get_CN_all(bond_range=3, printit=False)
        for cell in ([[0, 0, 0]] * 3, [[2, 0, 0], [4, 0, 0], [0, 0, 2]],
                     [[float("nan"), 0, 0], [0, 2, 0], [0, 0, 2]]):
            with self.subTest(cell=cell):
                cluster = make_cluster(Atoms("Cu", positions=[[0, 0, 0]], cell=cell, pbc=True))
                with self.assertRaises(ValueError):
                    cluster.get_CN_all(bond_range=3, printit=False)

    def test_empty_cluster_and_single_atom(self):
        empty = make_cluster(Atoms())
        self.assertEqual(empty.get_CN_all(printit=False), {})
        self.assertEqual(empty.get_CN_summary_all(printit=False), {})
        isolated = make_cluster(Atoms("Cu", positions=[[0, 0, 0]]))
        self.assertEqual(isolated.get_CN_for_atom_index(0, "Cu-Cu"), 0)
        self.assertEqual(isolated.get_CN_all(printit=False), {"Cu-Cu": {}})
        self.assertIsInstance(isolated.get_CN_summary_all(printit=False), dict)

    def test_valid_absent_species_returns_zero_and_invalid_symbol_rejects(self):
        cluster = make_cluster(Atoms("Cu2", positions=[[0, 0, 0], [2, 0, 0]]))
        self.assertEqual(cluster.get_CN("Cu", "O", printit=False)["Cu-O"], {})
        self.assertEqual(cluster.get_CN_for_atom_index(0, "Cu-O", cutoff=3), 0)
        self.assertEqual(cluster.get_CN_around_distance("Cu", "O", shell=1, printit=False), {0: 2})
        with self.assertRaises(ValueError):
            cluster.get_CN("Cu", "DefinitelyNotAnElement", printit=False)

    def test_atom_index_and_center_species_are_validated(self):
        cluster = make_cluster(Atoms("CuO", positions=[[0, 0, 0], [2, 0, 0]]))
        for index in (-1, 2, 0.5, True):
            with self.subTest(index=index), self.assertRaises((IndexError, TypeError, ValueError)):
                cluster.get_CN_for_atom_index(index, "Cu-O", cutoff=3)
        with self.assertRaises(ValueError):
            cluster.get_CN_for_atom_index(0, "O-Cu", cutoff=3)

    def test_natural_cutoff_uses_sum_of_atomic_radii(self):
        cluster = make_cluster(Atoms("Cu3", positions=[[0, 0, 0], [2.55, 0, 0], [20, 0, 0]]))
        self.assertEqual([cluster.get_CN_for_atom_index(i, "Cu-Cu") for i in range(3)], [1, 1, 0])

    def test_removal_includes_isolated_atoms_by_default_with_opt_out(self):
        cluster = make_cluster(Atoms("Cu3", positions=[[0, 0, 0], [2.55, 0, 0], [20, 0, 0]]))
        removed = quiet(cluster.remove_under_coordinated_atoms, "Cu", "Cu",
                        CN_threshold=2, bond_range=3)
        self.assertEqual(len(removed), 0)
        self.assertEqual(len(cluster.atoms), 3)
        retained = quiet(cluster.remove_under_coordinated_atoms, "Cu", "Cu", CN_threshold=2,
                         bond_range=3, include_isolated=False)
        np.testing.assert_allclose(retained.positions, [[20, 0, 0]])


class BackendParityTests(unittest.TestCase):
    def test_auto_falls_back_without_pymatgen(self):
        real_import = builtins.__import__

        def import_without_pymatgen(name, *args, **kwargs):
            if name.startswith("pymatgen"):
                raise ImportError("pymatgen intentionally unavailable for fallback test")
            return real_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=import_without_pymatgen):
            finite = make_cluster(Atoms("Cu2", positions=[[0, 0, 0], [2, 0, 0]]))
            self.assertEqual(finite.get_CN_all(bond_range=3, printit=False)["Cu-Cu"], {2: 1})
            self.assertEqual(finite.neighbor_backend_used, "scipy")
            periodic = make_cluster(bulk("Pt", "fcc", a=3.92))
            self.assertEqual(list(periodic.get_CN("Pt", "Pt", bond_range=3.2,
                                                 printit=False)["Pt-Pt"].values()), [12])
            self.assertEqual(periodic.neighbor_backend_used, "ase")

    def test_finite_backends_match_independent_atom_counts(self):
        rng = np.random.default_rng(1923)
        atoms = Atoms(["Cu", "O", "Pt"] * 10, positions=rng.uniform(-4, 4, (30, 3)))
        edges = [0, 1.5, 2.5, 4]
        backends = ["scipy", "ase", "auto"] + (["pymatgen"] if HAVE_PYMATGEN else [])
        for backend in backends:
            cluster = make_cluster(atoms.copy(), backend)
            cluster.get_CN_all(bond_range=4, shell_edges=edges, printit=False)
            for center in ("Cu", "O", "Pt"):
                for neighbor in ("Cu", "O", "Pt"):
                    with self.subTest(backend=backend, pair=f"{center}-{neighbor}"):
                        centers, expected = finite_shell_counts(atoms, center, neighbor, edges)
                        shells = cluster.CN_shells[f"{center}-{neighbor}"]
                        self.assertEqual(len(shells), 3)
                        for shell, counts in zip(shells, expected):
                            self.assertEqual(list(shell["atom_indices"]), centers)
                            np.testing.assert_array_equal(shell["atom_CN"], counts)
                            self.assertAlmostEqual(shell["average_CN"], np.mean(counts))
                            self.assertEqual(shell["distribution"],
                                             dict(zip(*np.unique(counts, return_counts=True))))

    @unittest.skipUnless(HAVE_PYMATGEN, "optional pymatgen backend is unavailable")
    def test_pymatgen_periodic_backend_matches_ase_multiple_images(self):
        atoms = bulk("Pt", "fcc", a=3.92)
        for backend in ("pymatgen", "ase"):
            with self.subTest(backend=backend):
                cluster = make_cluster(atoms.copy(), backend)
                cluster.get_CN_summary_all(bond_range=5, printit=False)
                self.assertEqual(list(cluster.CN["Pt-Pt"].values()), [12, 6, 24])
                self.assertEqual([shell["distribution"] for shell in cluster.CN_shells["Pt-Pt"]],
                                 [{12: 1}, {6: 1}, {24: 1}])


if __name__ == "__main__":
    unittest.main()
