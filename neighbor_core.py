import py3Dmol
import numpy as np
import plotly.graph_objects as go
from ase.io import read
from ase.io import write
import os
from scipy.spatial.distance import pdist
from ase.neighborlist import NeighborList, natural_cutoffs
from scipy.optimize import curve_fit
from collections import defaultdict
from ase.geometry import cell_to_cellpar
import pandas as pd
from itertools import product


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
    def __init__(self):
        """Initialize the ClusterNeighbor class.
        """
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
        self.cluster_size = 0
        self.indices = None

    def refresh_atoms(self):
        """ 
        When self.atoms is updated, refresh the parameters.
        """
        self._init_parameters()
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
        self.center = self.atoms.get_center_of_mass()
        
        # Create the XYZ string for visualization
        self.xyz_string = f"{len(self.atoms)}\n\n" 
        for atom in self.atoms:
            self.xyz_string += f"{atom.symbol} {atom.position[0]} {atom.position[1]} {atom.position[2]}\n"

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
        # Estimate the repetitions needed to ensure the desired diameter
        cell_lengths = self.atoms.get_cell_lengths_and_angles()[:3]
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
        if new_radius is None:
            new_radius = self.cluster_size - 0.1
        center = self.atoms[center_atom_index].position if center_atom_index is not None else self.center
        distances = np.linalg.norm(self.atoms.positions - center, axis=1)
        mask = distances <= new_radius
        atoms_smaller = self.atoms[mask]
        if self_apply:
            self.atoms = atoms_smaller
            self.refresh_atoms()
        return atoms_smaller

    def get_CN(self, center_atom=None, CN_atom=None, tolerance=0.01, bond_range=5, printit=True):
        """
        Calculate the coordination number (CN) for a specified bond type in the atomic cluster.
 
        Args:
        center_atom (str): Symbol of the central atom type for which CN is calculated. Defaults to the first atom's symbol.
        CN_atom (str): Symbol of the neighboring atom type for which CN is calculated. Defaults to the second atom's symbol.
        bond_ranges (float): Threshold for identifying significant gaps in bond lengths. Defaults to 0.01.
        tolerance (float): Maximum distance to consider for coordination number calculations. Defaults to 5.
        printit (bool): Whether to print the results. Defaults to True.

        Returns:
        dict: A dictionary where the keys are bond lengths and the values are the calculated coordination numbers.

        Example:
        """
        if center_atom is None:
            center_atom = list(self.pairs_types)[0].split('-')[0]
        if CN_atom is None:
            CN_atom = list(self.pairs_types)[0].split('-')[1]
            
        self.tolerance = tolerance
        self.bond_range = bond_range
        
        pairs_type = f"{center_atom}-{CN_atom}"
        center_atom_index = self.element_index_group[center_atom]
        CN_atom_index = self.element_index_group[CN_atom]

        # Create a NeighborList with the specified cutoff radius
        cutoffs = [bond_range] * len(self.atoms)
        nl = NeighborList(cutoffs, skin=0.5, bothways=True, self_interaction=False)
        nl.update(self.atoms)
        
        # Calculate all distances once and filter based on CN_atom_index
        distances_all = []
        for atom_i in center_atom_index:
            indices, offsets = nl.get_neighbors(atom_i)
            indices = indices[indices != atom_i]  # Exclude self-pairing
            if len(indices) == 0:
                continue
            distances = self.atoms.get_distances(atom_i, indices, mic=True)
            mask = np.isin(indices, CN_atom_index)
            distances_all.extend(distances[mask])

        # Convert the distances to a numpy array
        distances_all = np.array(distances_all)
        if len(distances_all) == 0:
            print(f"No valid distances found for bond type {pairs_type}.")
            return self.CN
    
        # Sort the distances and remove zeros and distances beyond the bond range
        distance_sorted = np.sort(distances_all)
        distance_sorted = distance_sorted[(distance_sorted != 0) & (distance_sorted < bond_range)]
        
        # Calculate differences and identify significant gaps
        if len(distance_sorted) == 0:
            print(f"No valid sorted distances found for bond type {pairs_type}.")
            return self.CN

        # Calculate differences and identify significant gaps
        diff = np.diff(distance_sorted)
        indices = np.where(diff > tolerance)[0] + 1
        self.CN_distances[pairs_type] = np.split(distance_sorted, indices)
        self.CN[pairs_type] = {np.average(group): group.shape[0] / self.elements_num[center_atom] for group in self.CN_distances[pairs_type]}
        
        if printit:
            print("=" * 20)
            print(pairs_type)
            print("=" * 20)
            for i, lengths_key in enumerate(self.CN[pairs_type].keys()):
                print(f"{i + 1} length: {lengths_key:.3f}  CN: {self.CN[pairs_type][lengths_key]:.3f}")

        return self.CN
    
    def get_CN_all(self, tolerance=0.01, bond_range=5, printit=True):
        """Calculate the coordination numbers for all pairs of atoms in the cluster.
        
        Args:
        tolerance (float): Threshold for identifying significant gaps in bond lengths. Defaults to 0.01.
        bond_range (float): Maximum distance to consider for coordination number calculations. Defaults to 5.
        printit (bool): Whether to print the results. Defaults to True.
        
        Returns:
        dict: A dictionary where the keys are bond types and the values are dictionaries containing the bond lengths and coordination numbers.
        """
        self.CN_distances = defaultdict(list)
        self.CN = {}
        for center_atom_i in self.elements_num.keys():
            for CN_atom_i in self.elements_num.keys():
                self.get_CN(center_atom=center_atom_i, 
                            CN_atom=CN_atom_i, 
                            tolerance=tolerance, 
                            bond_range=bond_range, 
                            printit=printit)
    
    def get_CN_around_distance(self, 
                               center_atom=None, 
                               CN_atom=None, 
                               target_distance=None, 
                               shell=None,
                               tolerance=0.01, 
                               bond_range=5, 
                               printit=True
                               ):
        """
        Calculate the coordination number (CN) for atoms with bond distances around the specified value.

        Args:
        target_distance (float): The target bond distance around which to count CNs. Defaults to None.
        pairs_type (str): The bond type for which CN is calculated (e.g., "Cu-Cu"). Defaults to None.
        tolerance (float): The tolerance for considering distances around the target distance. Defaults to 0.01.
        bond_range (float): Maximum distance to consider for coordination number calculations. Defaults to 5.
        printit (bool): Whether to print the results. Defaults to True.

        Returns:
        dict: A dictionary where the keys are coordination numbers and the values are the counts of atoms with that CN.
        """
        # If pairs_type is not provided, use the first pair type in the cluster
        if center_atom is None:
            center_atom = list(self.pairs_types)[0].split('-')[0]
        if CN_atom is None:
            CN_atom = list(self.pairs_types)[0].split('-')[1]
            
        pairs_type = f"{center_atom}-{CN_atom}"

        # If target_distance is not provided, use the CN data to get the target distance
        if target_distance is None:
            if self.CN.get(pairs_type) is None:
                self.get_CN(center_atom=center_atom, CN_atom=CN_atom, tolerance=tolerance, bond_range=bond_range, printit=printit)
            if shell is None:
                shell = 1
            target_distance = list(self.CN[pairs_type].keys())[shell-1]
        
        # Get the indices of the center and CN atoms
        center_atom_index = self.element_index_group[center_atom]
        CN_atom_index = self.element_index_group[CN_atom]

        # Create a NeighborList with the specified cutoff radius
        cutoffs = [bond_range] * len(self.atoms)
        nl = NeighborList(cutoffs, skin=tolerance, bothways=True, self_interaction=False)
        nl.update(self.atoms)
        
        # Count CNs for each atom based on distances around the target distance
        atom_CN_counts = {atom_i: 0 for atom_i in center_atom_index}
        for atom_i in center_atom_index:
            indices, offsets = nl.get_neighbors(atom_i)
            indices = indices[indices != atom_i]  # Exclude self-pairing
            if len(indices) == 0:
                continue
            distances = self.atoms.get_distances(atom_i, indices, mic=True)
            mask = (distances >= target_distance - tolerance) & (distances <= target_distance + tolerance) & np.isin(indices, CN_atom_index)
            atom_CN_counts[atom_i] = np.sum(mask)
        
        # Prepare the summary list
        CN_summary = defaultdict(int)
        for cn in atom_CN_counts.values():
            CN_summary[cn] += 1

        average_CN = np.mean(list(atom_CN_counts.values()))
        info = "\n".join([f"{count} atoms have CN of {cn}; " for cn, count in sorted(CN_summary.items())])

        self.CN_summary[pairs_type][target_distance] = {
            'average_CN': average_CN,
            'info': info,
            'tolerance': tolerance
        }

        if printit:
            print("=" * 20)
            print(f"Coordination numbers around {target_distance:.3f} Å for {pairs_type}")
            print("=" * 20)
            for cn, count in sorted(CN_summary.items()):
                print(f"{count} atoms have CN of {cn}")
                
        return CN_summary

    def get_CN_summary_all(self,
                        tolerance=0.01,
                        bond_range=5,
                        printit=True):
        """
        Calculate the coordination numbers for all pairs of atoms in the cluster.

        Args:
        tolerance (float): Threshold for identifying significant gaps in bond lengths. Defaults to 0.01.
        bond_range (float): Maximum distance to consider for coordination number calculations. Defaults to 5.
        printit (bool): Whether to print the results. Defaults to True.
        
        Returns:
        dict: A dictionary where the keys are bond types and the values are dictionaries containing the bond lengths and coordination numbers.
        """
        self.get_CN_all(tolerance=tolerance, bond_range=bond_range, printit=printit)
        for pair_type in self.pairs_types:
            center_atom, CN_atom = pair_type.split('-')
            for target_distance in self.CN[pair_type].keys():
                self.get_CN_around_distance(center_atom=center_atom,
                                            CN_atom=CN_atom,
                                            target_distance=target_distance,
                                            tolerance=tolerance,
                                            bond_range=bond_range,
                                            printit=printit)
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
        
    def remove_under_coordinated_atoms(self,
                                       center_atom=None,
                                       CN_atom=None, 
                                       CN_threshold=2, 
                                       bond_range=None, 
                                       self_apply=False):
        """ Remove under-coordinated atoms based on a specified coordination number threshold.

        Args:
            bond_type: The bond type for which the CN is evaluated (e.g., "Pb-O").
            CN_threshold: The threshold coordination number for identifying under-coordinated atoms.
            bond_range: The cutoff distance for considering neighbors.
            self_apply: Whether to apply the removal to the current cluster.

        Returns:
            atoms: The cluster with under-coordinated atoms removed.
        """
        if not hasattr(self, 'CN'):
            raise ValueError("CN data not available. Please run get_CN_all() first.")
        if center_atom is None:
            center_atom = list(self.pairs_types)[0].split('-')[0]
            print(f"Center atom not provided. Using {center_atom} as the center atom.")
        if CN_atom is None:
            CN_atom = list(self.pairs_types)[0].split('-')[1]
            print(f"CN atom not provided. Using {CN_atom} as the CN atom.")
            
        bond_type = f"{center_atom}-{CN_atom}"
        # Identify the indices of the atoms we are targeting for removal (atom1, e.g., "O" in "O-Pb")
        center_atom_indices = self.element_index_group[center_atom]
        
        # Collect indices of under-coordinated atoms
        under_coordinated_indices = []
        
        # Get the CN for each atom and remove under-coordinated atoms
        for atom_i in center_atom_indices:
            CN = self.get_CN_for_atom_index(atom_i, bond_type, cutoff=bond_range)     
            if CN > 0 and CN < CN_threshold:
                under_coordinated_indices.append(atom_i)
        
        if under_coordinated_indices:
            print(f"Remove under-coordinated atoms: {self.atoms[under_coordinated_indices]}[{under_coordinated_indices}]")
        else:
            print("No under-coordinated atoms found.")
            
        reduced_atoms = self.remove_atoms(under_coordinated_indices, self_apply=self_apply)
        
        return reduced_atoms
    
    def get_CN_for_atom_index(self, 
                        atom_index, 
                        bond_type, 
                        tolerance=0.01,
                        cutoff=None):
        """
        Extract the nearest neighbor coordination number (CN) for a specific atom index and bond type.

        Args:
        atom_index (int): The index of the atom for which the CN is calculated.
        bond_type (str): The bond type for which the CN is evaluated (e.g., "Pb-O").
        cutoff (float): The cutoff distance for considering neighbors. 

        Returns:
        int: The nearest neighbor coordination number for the specified atom and bond type.
        """
        # Extract atom types from the bond type
        atom1, atom2 = bond_type.split('-')

        # Set cutoffs if not provided
        cutoffs = natural_cutoffs(self.atoms) if cutoff is None else [cutoff] * len(self.atoms)
        if cutoff is None:
            cutoff = cutoffs[atom_index] 
        # Create a NeighborList with the specified cutoffs
        nl = NeighborList(cutoffs, skin=tolerance, bothways=True, self_interaction=False)
        nl.update(self.atoms)
        
        # Get neighbors and distances for the specified atom index
        indices, offsets = nl.get_neighbors(atom_index)
        distances = self.atoms.get_distances(atom_index, indices, mic=True)
        
        # Filter neighbors by bond type and find the nearest neighbors within the cutoff
        neighbor_symbols = np.array(self.atoms[indices].get_chemical_symbols())
        mask = (neighbor_symbols == atom2) & (distances <= cutoff)
        
        # Count the number of nearest neighbors matching the bond type
        cn = np.sum(mask)
        return int(cn)
    
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
        if not filename.endswith(".xlsx"):
            filename += ".xlsx"
        filepath = os.path.join(self.folder, filename)

        with pd.ExcelWriter(filepath, engine='openpyxl') as writer:
            for bond_type, distances in self.CN_summary.items():
                data = []
                for distance, details in distances.items():
                    N_atoms = details['info']
                    data.append({
                        "Bond length [A]": f"{distance:.3f}",
                        "range (tolerance) [A]": details['tolerance'],
                        "Average CN": details['average_CN'],
                        "N. atoms": N_atoms
                    })
                
                df = pd.DataFrame(data)
                df.to_excel(writer, sheet_name=bond_type, index=False)