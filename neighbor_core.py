import py3Dmol
import numpy as np
import plotly.graph_objects as go
from ase.io import read
from ase.io import write
import os
from scipy.spatial.distance import pdist
from ase.neighborlist import natural_cutoffs
from ase.data import atomic_numbers
from scipy.optimize import curve_fit
from collections import defaultdict
from ase.geometry import cell_to_cellpar
import pandas as pd
from itertools import product
from collections.abc import Mapping
import operator
from neighbor_search import BACKENDS, find_neighbor_edges


def fit_michaelis_menten(x, y, max_value=None, xlabel='x', ylabel='y', label_prefix="", maxfev=50000, extended_range=100):
    """
    Fit a Michaelis-Menten function with a constant offset. If max_value is provided, 
    the function is constrained to approach max_value as x -> +∞.
    """
    # Define the Michaelis-Menten function with a constant offset
    def michaelis_menten_function(x, alpha, beta, constant):
        return (alpha * x) / (1 + beta * x) + constant

    # Set initial guesses and bounds for parameters
    initial_guess = [np.max(y) / 2, 1, 0]  # Starting guess for alpha, beta, and constant

    if max_value is None:
        bounds = ([0, 0, -np.inf], [np.inf, np.inf, np.inf])  # Unbounded constant if max_value is not given
    else:
        assert max_value > 0, "max_value must be greater than zero"
        
        # Define a constrained fit function that adjusts constant based on max_value
        def constrained_function(x, alpha, beta):
            constant = max_value - (alpha / beta)
            return michaelis_menten_function(x, alpha, beta, constant)
        
        # Bounds only on alpha and beta; constant is calculated based on max_value
        bounds = ([0, 0], [np.inf, np.inf])
        initial_guess = [max_value / 2, 1]  # Adjusted initial guesses for alpha and beta

    # Perform the fitting
    if max_value is None:
        params, _ = curve_fit(michaelis_menten_function, x, y, p0=initial_guess, bounds=bounds, maxfev=maxfev)
        alpha, beta, constant = params
    else:
        params, _ = curve_fit(constrained_function, x, y, p0=initial_guess[:2], bounds=bounds, maxfev=maxfev)
        alpha, beta = params
        constant = max_value - (alpha / beta)

    # Create a DataFrame for the coefficients
    coeff_table = pd.DataFrame({
        'Coefficient': ['α', 'β', 'C'],
        'Value': [alpha, beta, constant]
    })

    # Generate extended and fitted data points
    x_extended = np.linspace(0, max(x) + extended_range, 100)
    y_extended = michaelis_menten_function(x_extended, alpha, beta, constant)
    x_fit = np.linspace(min(x), max(x), 100)
    y_fit = michaelis_menten_function(x_fit, alpha, beta, constant)

    # Create Plotly traces
    data_trace = go.Scatter(x=x, y=y, mode='markers', name=f'{label_prefix} Data Points')
    extended_fit_trace = go.Scatter(x=x_extended, y=y_extended, mode='lines', line=dict(dash='dash'), name=f'{label_prefix} Extended Fit')
    fit_trace = go.Scatter(x=x_fit, y=y_fit, mode='lines', name=f'{label_prefix} Michaelis-Menten Fit')

    return coeff_table, data_trace, extended_fit_trace, fit_trace

def fit_polynomial(x, y, degree, xlabel='x', ylabel='y', shell_label=""):
    """
    Fit a polynomial of the given degree to the provided points and return the plotly traces.
    """
    coeffs = np.polyfit(x, y, degree)
    polynomial = np.poly1d(coeffs)
    coeff_table = pd.DataFrame({
        'Coefficient': [f'x^{deg}' if deg > 0 else 'Constant' for deg in range(degree + 1)],
        'Value': coeffs[::-1]
    })

    x_fit = np.linspace(min(x), max(x), 100)
    y_fit = polynomial(x_fit)

    # Create traces instead of a complete figure
    data_trace = go.Scatter(x=x, y=y, mode='markers', name=f'{shell_label} Data Points')
    fit_trace = go.Scatter(x=x_fit, y=y_fit, mode='lines', name=f'{shell_label} Polynomial Fit (degree {degree})')

    return coeff_table, data_trace, fit_trace

class ClusterNeighbor:
    def __init__(self, neighbor_backend="auto", boundary_mode="auto"):
        """Analyze radial coordination shells without changing their definition.

        ``neighbor_backend`` is auto, pymatgen, scipy (finite only), or ase.
        Auto prefers pymatgen, falling back to SciPy/ASE when not installed.
        ``boundary_mode='auto'`` respects Atoms.pbc (including partial PBC);
        cluster disables PBC, and periodic enables all three directions.
        """
        self.neighbor_backend = neighbor_backend
        self.boundary_mode = boundary_mode
        self._validate_search_options()
        self.atoms = None
        self.folder = '.'
        self.xyz_string = ""
        self._init_parameters()

    def _init_parameters(self):
        """Initialize the parameters of the class.
        """
        self.center = None
        self.pairs_index = []
        self.pairs_element = []
        self.distance_all = []
        self.pairs_types = set()
        self.CN_distances = defaultdict(list)
        self.CN = {}
        self.CN_summary = defaultdict(lambda: defaultdict(dict))  # Initialize CN_summary as a nested defaultdict
        self.CN_shells = {}
        self._pair_parameters = {}
        self._neighbor_cache = None
        self.neighbor_backend_used = None
        self._atoms_signature = None
        self.pairs_group = {}
        self.cluster_size = 0
        self.indices = None

    def refresh_atoms(self):
        """ 
        When self.atoms is updated, refresh the parameters.
        """
        self._init_parameters()
        if self.atoms is None:
            raise ValueError("Load an ASE Atoms object before calculating neighbors")
        self.xyz_string = f"{len(self.atoms)}\n\n" + "\n".join(
            f"{atom.symbol} {atom.position[0]} {atom.position[1]} {atom.position[2]}" for atom in self.atoms
        )
        # get the elements and their counts            
        self.elements = self.atoms.get_chemical_symbols()
        unique_elements, counts = np.unique(self.elements, return_counts=True)
        self.elements_num = dict(zip(unique_elements, counts))
        
        # get pairs types
        self.pairs_types = set([f"{atom_i}-{atom_j}" for atom_i, atom_j in product(unique_elements, repeat=2)])

        # get the index of each element
        self.element_index_group = {element: [i for i, e in enumerate(self.elements) if e == element] for element in set(self.elements)}
        self.center = self.atoms.get_center_of_mass() if len(self.atoms) else np.zeros(3)
        self._symbols = np.asarray(self.elements)
        self._atoms_signature = self._structure_signature()

    def _validate_search_options(self):
        if self.neighbor_backend not in BACKENDS:
            raise ValueError("neighbor_backend must be auto, pymatgen, scipy or ase")
        if self.boundary_mode not in {"auto", "cluster", "periodic"}:
            raise ValueError("boundary_mode must be auto, cluster or periodic")

    def _structure_signature(self):
        return (self.atoms.positions.tobytes(), self.atoms.numbers.tobytes(),
                self.atoms.cell.array.tobytes(), self.atoms.pbc.tobytes(),
                self.atoms.get_masses().tobytes(), self.neighbor_backend,
                self.boundary_mode)

    def _sync_atoms(self):
        """Invalidate derived results on the next calculation after direct edits."""
        self._validate_search_options()
        if self.atoms is None:
            raise ValueError("Load an ASE Atoms object before calculating neighbors")
        if self._atoms_signature != self._structure_signature():
            self.refresh_atoms()

    def _get_neighbor_edges(self, radius):
        # Public callers synchronize once before requesting one or many pairs.
        if self._neighbor_cache is None or self._neighbor_cache[0] < radius:
            pbc = self.atoms.pbc
            if self.boundary_mode != "auto":
                pbc = np.full(3, self.boundary_mode == "periodic")
            result = find_neighbor_edges(self.atoms, radius, pbc, self.neighbor_backend)
            self._neighbor_cache = (radius, result)
            self.neighbor_backend_used = result[-1]
        return self._neighbor_cache[1][:4]

    def _calculate_distance(self, coord1, coord2):
        """Calculate the Euclidean distance between two sets of coordinates.

        Args:
        coord1 (list): The first set of coordinates.
        coord2 (list): The second set of coordinates.
        
        Returns:
        float: The Euclidean distance between the two sets of coordinates.
        """
        return np.linalg.norm(np.array(coord1) - np.array(coord2))

    def load_xyz(self, from_file=True, path=None, atom_object=None):
        """Load the atomic cluster from an XYZ file or an ASE Atoms object.

        Args:
        from_file (bool): Whether to load the cluster from a file. Defaults to True.
        path (str): The path to the XYZ file. Defaults to None.
        atom_object (Atoms): The ASE Atoms object representing the atomic cluster. Defaults to None.
        
        Raises:
        ValueError: If neither a file path nor an ASE Atoms object is provided.
        
        Returns:
        Atoms: The ASE Atoms object representing the atomic cluster.
        """
        # Load the atoms from a file or an object
        if from_file:
            self.atoms = read(path)
            self.folder = os.path.dirname(path)
        else:
            self.atoms = atom_object
        
        self.refresh_atoms()
    
    def save(self, filename):
        """Save the atomic cluster to an XYZ file.

        Args:
        filename (str): The name of the file to save the atomic cluster to.
        
        Returns:
        Atoms: The ASE Atoms object representing the atomic cluster.
        """
        write(filename, self.atoms)

    def expand_cif(self, replication_factors=(2, 2, 2), self_apply=False):
        """ Expand the cluster by replicating it in all three dimensions.

        Args:
            replication_factors (tuple): The number of times to replicate the cluster in each dimension.
            self_apply (bool): Whether to apply the expansion to the current cluster.
        Returns:
            expanded_cluster (Atoms): The expanded cluster.
        """
        expanded_cluster = self.atoms.repeat(replication_factors)
        if self_apply:
            self.atoms = expanded_cluster 
            self.refresh_atoms()
        return expanded_cluster


    def expand_to_sphere(self, target_diameter=50, self_apply=False):
        """Expand the CIF structure to cover a desired diameter and cut a spherical cluster around the center of mass.

        Args:
            target_diameter (float): The desired size of the cluster in Ångströms.
            self_apply (bool): Whether to apply the expansion and cut to the current cluster.

        Returns:
            expanded_cluster (Atoms): The expanded and cut cluster.
        """
        self._sync_atoms()
        _, target_diameter = self._validate_cn_parameters(0., target_diameter)
        # Estimate the repetitions needed to ensure the desired diameter
        cell_lengths = self.atoms.cell.lengths()
        if np.any(cell_lengths <= 0):
            raise ValueError("Expanding a crystal requires three nonzero cell vectors")
        reps = np.ceil(target_diameter / np.array(cell_lengths)).astype(int)

        # Expand the unit cell using atoms.repeat to cover the desired diameter
        expanded_atoms = self.atoms.repeat(reps)

        # Recalculate the center of mass for the expanded structure
        center_of_mass_expanded = expanded_atoms.get_center_of_mass()

        # Calculate the distance of each atom from the center of mass
        distances = np.linalg.norm(expanded_atoms.positions - center_of_mass_expanded, axis=1)

        # Select atoms within the desired spherical diameter
        radius = target_diameter / 2.0
        mask = distances <= radius
        expanded_cluster = expanded_atoms[mask]
        expanded_cluster.set_pbc(False)

        # Sort the atoms by atomic number
        sorted_indices = np.argsort(expanded_cluster.get_atomic_numbers())
        expanded_cluster = expanded_cluster[sorted_indices]

        # Apply the expansion and cut to the current cluster if self_apply is True
        if self_apply:
            self.atoms = expanded_cluster
            self.refresh_atoms()

        return expanded_cluster


    def view_xyz(self, style_all=None, highlight_atom1="O", highlight_atom2="Pb", label=False, show_symbol=False):
        """Visualize the atomic cluster using py3Dmol.

        Args:
            style_all (_type_, optional): _description_. Defaults to None.
            highlight_atom1 (str, optional): _description_. Defaults to "O".
            highlight_atom2 (str, optional): _description_. Defaults to "Pb".
            label (bool, optional): _description_. Defaults to False.
            show_symbol (bool, optional): _description_. Defaults to False.
        """
        view = py3Dmol.view(width=400, height=400)
        
        # Add the model to the view
        view.addModel(self.xyz_string, 'xyz')
        
        if style_all is None:
            style_all = {'stick': {'radius': .1, 'alpha': 0.2, 'color': 'gray'}, 
                         'sphere': {'radius': .3}}
        view.setStyle(style_all)
        view.addStyle({'atom': highlight_atom1}, {'sphere': {'color': 'red', 'radius': 0.5}})
        view.addStyle({'atom': highlight_atom2}, {'sphere': {'color': 'blue', 'radius': 0.3}})
        view.setBackgroundColor('0xeeeeee')

        if label:
            for i, atom in enumerate(self.atoms):
                symbol = atom.symbol if show_symbol else ""
                view.addLabel(f"{i}{symbol}", {'position': {'x': atom.position[0], 'y': atom.position[1], 'z': atom.position[2]}, 
                                    'fontColor': 'k', 'fontSize': 12, 'backgroundColor': 'white', 'backgroundOpacity': 0.5})
        
        view.zoomTo()
        view.show()
        view.title(self.atoms.get_chemical_formula())

        # Assign the new view to the instance variable
        self.view = view

    def get_cluster_size(self):
        """Estimate the cluster size using the maximum pairwise distance method.

        Returns:
        float: Approximate cluster size in diameter.
        """
        self._sync_atoms()
        if len(self.atoms) < 2:
            self.cluster_size = 0.
            return self.cluster_size
        positions = self.atoms.get_positions()
        distances = pdist(positions)
        self.cluster_size = np.max(distances)
        return self.cluster_size
    
    def get_cluster_size_bounding_box(self):
            """
            Estimate the cluster size using the bounding box method.

            Returns:
            float: Approximate cluster size.
            """
            cellpar = cell_to_cellpar(self.atoms.get_cell())
            max_dimension = max(cellpar[:3])
            return max_dimension

    def shrink_cluster_size(self, new_radius=None, center_atom_index=None, self_apply=False):
        """Shrink the cluster size by removing atoms outside a specified radius.

        Args:
            new_radius (_type_, optional): _description_. Defaults to None.
            center_atom_index (_type_, optional): _description_. Defaults to None.
            self_apply (bool, optional): _description_. Defaults to False.

        Returns:
            atoms: The shrunken cluster.
        """
        self._sync_atoms()
        if new_radius is None:
            new_radius = self.cluster_size - 0.1
        center = self.atoms[center_atom_index].position if center_atom_index is not None else self.center
        distances = np.linalg.norm(self.atoms.positions - center, axis=1)
        mask = distances <= new_radius
        atoms_smaller = self.atoms[mask]
        atoms_smaller.set_pbc(False)
        if self_apply:
            self.atoms = atoms_smaller
            self.refresh_atoms()
        return atoms_smaller

    @staticmethod
    def _validate_cn_parameters(tolerance, bond_range):
        tolerance, bond_range = float(tolerance), float(bond_range)
        if not np.isfinite(tolerance) or tolerance < 0:
            raise ValueError("tolerance must be finite and nonnegative")
        if not np.isfinite(bond_range) or bond_range <= 0:
            raise ValueError("bond_range must be finite and positive")
        return tolerance, bond_range

    def _resolve_pair(self, center_atom, CN_atom):
        if center_atom is None:
            if not self.elements:
                raise ValueError("Specify elements when querying an empty cluster")
            center_atom = self.elements[0]
        if CN_atom is None:
            CN_atom = self.elements[min(1, len(self.elements) - 1)] if self.elements else center_atom
        for symbol in (center_atom, CN_atom):
            if symbol not in atomic_numbers:
                raise ValueError(f"Unknown element symbol: {symbol!r}")
        return center_atom, CN_atom, f"{center_atom}-{CN_atom}"

    @staticmethod
    def _resolve_shell_edges(shell_edges, center_atom, CN_atom, bond_range):
        if isinstance(shell_edges, Mapping):
            # A single A-B definition also covers B-A unless explicitly overridden.
            keys = (f"{center_atom}-{CN_atom}", (center_atom, CN_atom),
                    f"{CN_atom}-{center_atom}", (CN_atom, center_atom))
            shell_edges = next((shell_edges[k] for k in keys if k in shell_edges), None)
        if shell_edges is None:
            return None
        edges = np.asarray(shell_edges, dtype=float)
        if (edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all()
                or edges[0] < 0 or np.any(np.diff(edges) <= 0)):
            raise ValueError("shell_edges must be a finite, nonnegative, strictly increasing sequence")
        if edges[-1] > bond_range:
            raise ValueError("The final shell edge must not exceed bond_range")
        return tuple(float(value) for value in edges)

    def _pair_edges(self, center_atom, CN_atom, bond_range):
        i, j, _, distances = self._get_neighbor_edges(bond_range)
        mask = ((self._symbols[i] == center_atom) & (self._symbols[j] == CN_atom)
                & (distances < bond_range))
        return i[mask], distances[mask]

    def _shell_record(self, center_atom, centers, distances, shell_id, lower, upper, method):
        atom_indices = np.asarray(self.element_index_group.get(center_atom, []), dtype=int)
        counts = np.bincount(centers, minlength=len(self.atoms))[atom_indices]
        values, frequencies = np.unique(counts, return_counts=True)
        distribution = dict(zip(map(int, values), map(int, frequencies)))
        return {
            'shell_id': shell_id,
            'distance': float(np.mean(distances)) if len(distances) else None,
            'lower_bound': lower, 'upper_bound': upper, 'method': method,
            'distance_min': float(np.min(distances)) if len(distances) else None,
            'distance_max': float(np.max(distances)) if len(distances) else None,
            'distance_std': float(np.std(distances)) if len(distances) else None,
            # An absent center species has no defined mean, rather than a zero CN.
            'average_CN': float(np.mean(counts)) if len(counts) else None,
            'atom_indices': atom_indices, 'atom_CN': counts,
            'distribution': distribution, 'distances': distances,
        }

    def _ensure_pair_shells(self, center_atom, CN_atom, tolerance, bond_range, shell_edges):
        pair = f"{center_atom}-{CN_atom}"
        edges = self._resolve_shell_edges(shell_edges, center_atom, CN_atom, bond_range)
        parameters = (tolerance, bond_range, edges)
        if self._pair_parameters.get(pair) == parameters:
            return self.CN_shells[pair]
        centers, distances = self._pair_edges(center_atom, CN_atom, bond_range)
        order = np.argsort(distances)
        centers, distances = centers[order], distances[order]
        records = []
        if edges is None:
            if len(distances):
                bounds = np.r_[0, np.flatnonzero(np.diff(distances) > tolerance) + 1, len(distances)]
                for shell_id, (lo, hi) in enumerate(zip(bounds[:-1], bounds[1:]), 1):
                    records.append(self._shell_record(
                        center_atom, centers[lo:hi], distances[lo:hi], shell_id,
                        float(distances[lo]), float(distances[hi - 1]), 'gap'))
        else:
            for shell_id, (lo, hi) in enumerate(zip(edges[:-1], edges[1:]), 1):
                start, stop = np.searchsorted(distances, [lo, hi], side='left')
                records.append(self._shell_record(
                    center_atom, centers[start:stop], distances[start:stop], shell_id, lo, hi, 'edges'))
        self.CN_shells[pair] = records
        self.CN_distances[pair] = [record['distances'] for record in records]
        self.CN[pair] = {record['distance']: record['average_CN']
                         for record in records if record['distance'] is not None}
        self.CN_summary[pair] = {}
        self._pair_parameters[pair] = parameters
        return records

    @staticmethod
    def _summary_record(record, tolerance):
        result = {key: value for key, value in record.items() if key != 'distances'}
        result['tolerance'] = tolerance
        result['info'] = "\n".join(
            f"{count} atoms have CN of {cn}; " for cn, count in record['distribution'].items())
        return result

    def get_CN(self, center_atom=None, CN_atom=None, tolerance=0.01, bond_range=5,
               printit=True, *, shell_edges=None):
        """Calculate mean radial-shell CN for one directed element pair.

        All central atoms, including those with zero neighbors, enter the mean.
        Distances are strictly less than ``bond_range`` (Angstrom).
        By default, sorted distances are split where adjacent gaps exceed
        ``tolerance``. For distorted structures or size series, ``shell_edges``
        sets fixed [lower, upper) intervals. It may be a sequence or a mapping
        from element pair to sequence; a reverse-pair definition is reused.
        Unspecified pairs use gap grouping. Empty explicit shells are retained
        in ``CN_shells``; the legacy distance-keyed ``CN`` contains occupied ones.
        """
        self._sync_atoms()
        tolerance, bond_range = self._validate_cn_parameters(tolerance, bond_range)
        center_atom, CN_atom, pair = self._resolve_pair(center_atom, CN_atom)
        self.tolerance, self.bond_range = tolerance, bond_range
        self._ensure_pair_shells(center_atom, CN_atom, tolerance, bond_range, shell_edges)
        if printit:
            print("=" * 20, pair, "=" * 20, sep="\n")
            for record in self.CN_shells[pair]:
                if record['distance'] is not None:
                    print(f"{record['shell_id']} length: {record['distance']:.3f}  CN: {record['average_CN']:.3f}")
        return self.CN

    def get_CN_all(self, tolerance=0.01, bond_range=5, printit=True, *, shell_edges=None):
        """Calculate every directed element pair using one shared neighbor table."""
        self._sync_atoms()
        tolerance, bond_range = self._validate_cn_parameters(tolerance, bond_range)
        self.CN_distances = defaultdict(list)
        self.CN = {}
        self.CN_shells = {}
        self.CN_summary = defaultdict(lambda: defaultdict(dict))
        self._pair_parameters = {}
        for center_atom in self.elements_num:
            for CN_atom in self.elements_num:
                self.get_CN(center_atom, CN_atom, tolerance, bond_range, printit,
                            shell_edges=shell_edges)
        return self.CN

    def get_CN_around_distance(self, center_atom=None, CN_atom=None,
                               target_distance=None, shell=None, tolerance=0.01,
                               bond_range=5, printit=True, *, shell_edges=None):
        """Return a histogram mapping per-atom CN to number of central atoms.

        ``shell`` is 1-based and counts the exact members used by ``get_CN``;
        the default is shell 1. Supplying ``target_distance`` instead explicitly
        requests the inclusive window target +/- tolerance, restricted to
        distances < bond_range. The two modes are intentionally distinct.
        """
        self._sync_atoms()
        tolerance, bond_range = self._validate_cn_parameters(tolerance, bond_range)
        center_atom, CN_atom, pair = self._resolve_pair(center_atom, CN_atom)
        if target_distance is not None and shell is not None:
            raise ValueError("Specify either shell or target_distance, not both")
        records = self._ensure_pair_shells(center_atom, CN_atom, tolerance, bond_range, shell_edges)
        if target_distance is None:
            shell = 1 if shell is None else operator.index(shell)
            if shell < 1 or (records and shell > len(records)) or (not records and shell != 1):
                raise ValueError(f"Shell {shell} is unavailable for {pair}")
            if not records:
                count = len(self.element_index_group.get(center_atom, []))
                return {0: count} if count else {}
            record = records[shell - 1]
            key = record['distance']
        else:
            target_distance = float(target_distance)
            if not np.isfinite(target_distance) or target_distance < 0:
                raise ValueError("target_distance must be finite and nonnegative")
            centers, distances = self._pair_edges(center_atom, CN_atom, bond_range)
            lower, upper = max(0., target_distance - tolerance), target_distance + tolerance
            mask = (distances >= lower) & (distances <= upper)
            record = self._shell_record(center_atom, centers[mask], distances[mask],
                                        None, lower, min(upper, bond_range), 'window')
            key = target_distance
        if key is not None:
            self.CN_summary[pair][key] = self._summary_record(record, tolerance)
        if printit:
            label = f"shell {record['shell_id']}" if target_distance is None else f"window around {target_distance:.3f} Å"
            print(f"Coordination numbers for {pair}, {label}")
            for cn, count in record['distribution'].items():
                print(f"{count} atoms have CN of {cn}")
        return dict(record['distribution'])

    def get_CN_summary_all(self, tolerance=0.01, bond_range=5, printit=True, *, shell_edges=None):
        """Summarize the same shell memberships used for mean CN calculation.

        Legacy distance-keyed summaries omit empty explicit shells. Full records,
        including empty intervals and stable shell IDs, are in ``CN_shells``.
        """
        self.get_CN_all(tolerance, bond_range, printit=False, shell_edges=shell_edges)
        for pair, records in self.CN_shells.items():
            for record in records:
                if record['distance'] is not None:
                    self.CN_summary[pair][record['distance']] = self._summary_record(record, float(tolerance))
        if printit:
            self.print_CN_summary()
        return self.CN_summary

    def get_pairs(self):
        """Get all pairs of atoms in the cluster and their distances.

        Returns:
        dict: A dictionary containing the pairs of atoms and their distances.
        """
        num_atoms = len(self.atoms)
        pairs_index = np.array(np.triu_indices(num_atoms, k=1)).T
        
        # Get the symbols of the atoms in each pair        
        symbols = np.array(self.atoms.get_chemical_symbols())
        pairs_element = symbols[pairs_index]

        pairs = np.array([f"{atom_i}-{atom_j}" for atom_i, atom_j in pairs_element])

        distance_matrix = self.atoms.get_all_distances()
        distance_all = distance_matrix[pairs_index[:, 0], pairs_index[:, 1]]

        self.pairs_index = pairs_index.tolist()
        self.pairs_element = pairs_element.tolist()
        self.pairs = pairs.tolist()
        self.distance_all = distance_all.tolist()

        self.pairs_group = {
            key: {
                'pairs_index': pairs_index[np.array(self.pairs) == key].tolist(),
                'distance': distance_all[np.array(self.pairs) == key].tolist()
            } for key in self.pairs_types
        }

        return self.pairs_group
    
    def plot_hist(self, binsize=0.2, plot_engine="plt"):
        """Plot a histogram of the distances between pairs of atoms in the cluster.

        Args:
        binsize (float): The size of the bins for the histogram. Defaults to 0.2.
        plot_engine (str): The plotting engine to use ('plt' or 'plotly'). Defaults to 'plt'.
        """
        if len(self.pairs_group) == 0:
            self.get_pairs()
        
        if plot_engine == "plotly":
            fig = go.Figure()
            for key_i in self.pairs_group.keys():
                fig.add_trace(go.Histogram(x=self.pairs_group[key_i]['distance'], name=key_i, opacity=0.6, 
                                        xbins={'size':binsize},marker={'line':{'color':'white','width':2}}))

            fig.update_layout(
                xaxis_title_text='Distances [A]', yaxis_title_text='pairs',
                plot_bgcolor='rgba(0.02,0.02,0.02,0.02)',  # Transparent plot background
                xaxis={'tickmode':'auto'}, barmode='overlay',  # Overlay histograms,
                width=600, height=400)
            fig.show()
            
        elif plot_engine == "plt":
            plt.figure(figsize=(8,3))
            for key_i in self.pairs_group.keys():
                plt.hist(self.pairs_group[key_i]['distance'], bins=80, alpha=0.3, edgecolor='white', label=key_i)
            plt.xlabel("Distance [Å]")
            plt.ylabel("Number of pairs")
            plt.legend()
            plt.tight_layout()
            plt.show()
    
    def remove_atoms(self, indices_remove_lst=[-1], self_apply=False):
        indices_keep_lst = [i for i in range(len(self.atoms)) if i not in indices_remove_lst]
        reduced_atoms = self.atoms[indices_keep_lst]
        if self_apply:
            self.atoms = reduced_atoms
            self.refresh_atoms()
        return reduced_atoms
        
    def _coordination_counts(self, center_atom, CN_atom, cutoff):
        """Counts for all atom indices, using a scalar or summed natural radii."""
        if cutoff is None:
            radii = np.asarray(natural_cutoffs(self.atoms), dtype=float)
            centers = self.element_index_group.get(center_atom, [])
            neighbors = self.element_index_group.get(CN_atom, [])
            if not centers or not neighbors:
                return np.zeros(len(self.atoms), dtype=int)
            radius = float(radii[centers].max() + radii[neighbors].max())
            if radius <= 0:
                return np.zeros(len(self.atoms), dtype=int)
        else:
            _, radius = self._validate_cn_parameters(0., cutoff)
        i, j, _, distances = self._get_neighbor_edges(radius)
        limits = radii[i] + radii[j] if cutoff is None else radius
        mask = ((self._symbols[i] == center_atom) & (self._symbols[j] == CN_atom)
                & (distances <= limits))
        return np.bincount(i[mask], minlength=len(self.atoms))

    def remove_under_coordinated_atoms(self, center_atom=None, CN_atom=None,
                                       CN_threshold=2, bond_range=None,
                                       self_apply=False, *, include_isolated=True):
        """Remove atoms with CN below the threshold in one batch.

        Includes CN=0 by default; use ``include_isolated=False`` to preserve
        isolated atoms. This is a single pass on the original structure, not
        iterative pruning. ``bond_range=None`` uses summed covalent radii.
        """
        self._sync_atoms()
        center_atom, CN_atom, _ = self._resolve_pair(center_atom, CN_atom)
        threshold = float(CN_threshold)
        if not np.isfinite(threshold) or threshold < 0:
            raise ValueError("CN_threshold must be finite and nonnegative")
        counts = self._coordination_counts(center_atom, CN_atom, bond_range)
        centers = np.asarray(self.element_index_group.get(center_atom, []), dtype=int)
        selected = counts[centers] < threshold
        if not include_isolated:
            selected &= counts[centers] > 0
        indices = centers[selected].tolist()
        if indices:
            print(f"Remove under-coordinated atoms: {self.atoms[indices]}[{indices}]")
        else:
            print("No under-coordinated atoms found.")
        return self.remove_atoms(indices, self_apply=self_apply)

    def get_CN_for_atom_index(self, atom_index, bond_type, tolerance=0.01, cutoff=None):
        """Return integer CN for an atom and directed element pair.

        A supplied ``cutoff`` is an inclusive pair distance. Otherwise use the
        sum of the two atoms' covalent radii. ``tolerance`` is retained for API
        compatibility and does not broaden that scientific distance threshold.
        Atom indices must be nonnegative integers and match the central element.
        """
        self._sync_atoms()
        if isinstance(atom_index, (bool, np.bool_)):
            raise ValueError("atom_index must be a nonnegative integer")
        atom_index = operator.index(atom_index)
        if atom_index < 0 or atom_index >= len(self.atoms):
            raise IndexError("atom_index is outside the cluster")
        self._validate_cn_parameters(tolerance, 1.)
        if not isinstance(bond_type, str) or len(bond_type.split('-')) != 2:
            raise ValueError("bond_type must be a directed pair such as 'Cu-O'")
        center_atom, CN_atom, _ = self._resolve_pair(*bond_type.split('-'))
        if self.elements[atom_index] != center_atom:
            raise ValueError("The atom's element does not match the central element in bond_type")
        counts = self._coordination_counts(center_atom, CN_atom, cutoff)
        return int(counts[atom_index])

    def print_CN_summary(self):
        for pairs_type, distances in self.CN_summary.items():
            print("=" * 35)
            print(f"Bond Type: {pairs_type}")
            print("=" * 35)
            for distance, details in distances.items():
                print(f"-- Distance: {distance}")
                print(f"-- Average CN: {details['average_CN']}")
                print(f"-- Info: \n{details['info']}")
                print(f"-- Tolerance: {details['tolerance']} \n")
    
    def write_to_excel(self, filename="output"):
        """Write the coordination number data to an Excel file.

        Args:
        filename (str): The name of the Excel file to write the data to. Defaults to "output.xlsx".
        """
        self._sync_atoms()
        filename = os.fspath(filename)
        if not filename.endswith(".xlsx"):
            filename += ".xlsx"
        filepath = os.path.join(self.folder, filename)

        with pd.ExcelWriter(filepath, engine='openpyxl') as writer:
            pairs = sorted(set(self.CN_shells) | set(self.CN_summary))
            for bond_type in pairs:
                distances = self.CN_summary.get(bond_type, {})
                records = self.CN_shells.get(bond_type, [])
                tolerance = self._pair_parameters.get(bond_type, (None,))[0]
                # Retain explicit window queries, and include empty fixed shells
                # so exported size series keep their original shell IDs.
                rows = list(distances.items())
                if not rows:
                    rows = [(record['distance'], self._summary_record(record, tolerance))
                            for record in records]
                else:
                    rows.extend((None, self._summary_record(record, tolerance))
                                for record in records if record['distance'] is None)
                rows.sort(key=lambda row: (row[1].get('shell_id') is None,
                                          row[1].get('shell_id') or 0,
                                          row[0] if row[0] is not None else -1))
                data = []
                for distance, details in rows:
                    N_atoms = details['info']
                    data.append({
                        "Bond length [A]": f"{distance:.3f}" if distance is not None else None,
                        "range (tolerance) [A]": details['tolerance'],
                        "Average CN": details['average_CN'],
                        "N. atoms": N_atoms,
                        "Shell": details.get('shell_id'),
                        "Shell method": details.get('method'),
                        "Lower bound [A]": details.get('lower_bound'),
                        "Upper bound [A]": details.get('upper_bound'),
                        "Distance std [A]": details.get('distance_std'),
                    })
                
                df = pd.DataFrame(data)
                df.to_excel(writer, sheet_name=bond_type, index=False)
            if not pairs:
                pd.DataFrame(columns=["Bond length [A]", "Average CN"]).to_excel(
                    writer, sheet_name="Summary", index=False)
