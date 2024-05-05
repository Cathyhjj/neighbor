
# -----------------------------------------------------------------------------
# :author:    Juanjuan Huang
# :email:     juanjuan.huang@anl.gov
# :copyright: Copyright © 2023, UChicago Argonne, LLC
# -----------------------------------------------------------------------------

import py3Dmol
import numpy as np
import matplotlib.pyplot as plt
# from IPython.display import display
from ipyfilechooser import FileChooser
# import re
import ase
from ase.io import read
from ase.io import write
import plotly.graph_objects as go
import copy
from pprint import pprint
import ipywidgets as widgets
from IPython.display import display, clear_output
import os

class ClusterNeighbor(object):
    def __init__(self):
        pass
    
    def _calculate_distance(self, coord1, coord2):
        return np.sqrt((coord2[0] - coord1[0])**2 + (coord2[1] - coord1[1])**2 + (coord2[2] - coord1[2])**2)


    def load_xyz(self, from_file=True, path=None,  atom_object=None):
        if from_file is True:
            self.atoms = read(path)
        else:
            self.atoms = atom_object
        self.elements = [self.atoms[i].symbol for i in range(len(self.atoms))]
        self.elements_num = {element_i: self.elements.count(element_i) for element_i in set(self.elements)}
        self.element_index_group = {element_set_i: [i for i, element_i in enumerate(self.elements) if element_i == element_set_i] for element_set_i in set(self.elements)}
        self.center = self.atoms.get_center_of_mass()
        
        self.xyz_string = f"{len(self.atoms)}\n\n" 
        for atom in self.atoms:
            self.xyz_string += f"{atom.symbol} {atom.position[0]} {atom.position[1]} {atom.position[2]}\n"
        
    def expand_cif(self, replication_factors = (2, 2, 2)):
        expanded_cluster = self.atoms.repeat(replication_factors)
        return expanded_cluster
    
    def view_xyz(self, style_all=None, highlight_atom1="O", highlight_atom2="Pb", label=True):
        self.view = py3Dmol.view(width=500,height=500)
        self.view.addModel(self.xyz_string,'xyz',)
        if style_all is None:
            style_all = {'stick':{'radius':.1, 'alpha':0.2, 'color':'gray'}, 
                         'sphere': {'radius':.3}
                        }
        self.view.setStyle(style_all)
        
        self.view.addStyle({'atom': highlight_atom1}, 
                           {'sphere': {'color': 'red', 'radius': 0.5}})  
        
        self.view.addStyle({'atom': highlight_atom2}, 
                           {'sphere': {'color': 'blue', 'radius': 0.3}})  
        self.view.setBackgroundColor('0xeeeeee')
        if label:
            for i, atom_i in enumerate(self.atoms):
                self.view.addLabel(f"{i}", {'position': {'x': atom_i.position[0], 'y': atom_i.position[1], 'z': atom_i.position[2]}, 
                                    'fontColor': 'k', 'fontSize': 12, 'backgroundColor': 'white', 'backgroundOpacity':0.5})
        self.view.zoomTo()
        self.view.show()
        self.view.title(self.atoms.get_chemical_formula())
    
    def get_cluster_size(self):
        self.cluster_size = self.atoms.get_all_distances().max()/2
        # print(f"Cluster size is {self.cluster_size} A")
        return self.cluster_size
    
    def shrink_cluster_size(self, new_radius=None):
        if new_radius is None:
            new_radius = self.cluster_size - 0.1
        
        atoms_smaller = copy.deepcopy(self.atoms)
        indices_remove_lst = []
        for i, atom_i in enumerate(atoms_smaller):
            radius_i = np.abs(self._calculate_distance(self.atoms.get_positions()[i], self.center))
            if radius_i > new_radius:
                indices_remove_lst.append(i)

        atoms_smaller = self.remove_atoms(indices_remove_lst)
        return atoms_smaller
    
    def remove_atoms(self, indices_remove_lst=[-1]):
        indices_keep_lst = [i for i in range(len(self.atoms)) if i not in indices_remove_lst]
        atoms_new = self.atoms[indices_keep_lst]
        return atoms_new

    def get_pairs(self):
        self.pairs_index = [(i, j) for i in range(len(self.atoms)) for j in range(i + 1, len(self.atoms))]
        # self.pairs_index = [(i, j) for i in range(len(self.atoms)) for j in range(len(self.atoms))]
        self.pairs_element = [sorted([self.atoms[i].symbol, self.atoms[j].symbol]) for i, j in self.pairs_index]
        self.pairs = [f"{atom_i}-{atom_j}" for atom_i, atom_j in self.pairs_element]
        self.pairs_unique = [f"{self.atoms[i].symbol}({self.atoms[i].index})-{self.atoms[j].symbol}({self.atoms[j].index})" for i, j in self.pairs_index]
        self.distance_all = [self.atoms.get_all_distances()[i][j] for i, j in self.pairs_index]
        self.pairs_types = set(self.pairs)
        
        self.pairs_group = {key: {'pairs_index':[],
                                  'pairs':[],
                                  'pairs_unique':[],
                                  'distance':[]} for key in self.pairs_types}
        
        for i in range(len(self.pairs)):
            self.pairs_group[self.pairs[i]]['pairs_index'].append(self.pairs_index[i])
            self.pairs_group[self.pairs[i]]['pairs_unique'].append(self.pairs_unique[i])
            self.pairs_group[self.pairs[i]]['distance'].append(self.distance_all[i])
        return self.pairs_group

    def get_CN(self, center_atom=None, CN_atom=None, error_bar=0.01, CN_range=5, printit=True):
        if not hasattr(self, 'pairs_types'):
            self.get_pairs()
            
        # if not hasattr(self, 'CN_distances'):
        self.CN_distances = {}
        self.CN = {}
                    
        if center_atom is None:
            center_atom = self.atoms[0].symbol
        
        if CN_atom is None:
            CN_atom = self.atoms[1].symbol
            
        bond_type = f"{center_atom}-{CN_atom}"                
        center_atom_index = self.element_index_group[center_atom]
        CN_atom_index = self.element_index_group[CN_atom]
        
        distances_all = np.asarray([self.atoms.get_distances(atom_i, CN_atom_index) for atom_i in center_atom_index])
        distance_sorted = np.sort(distances_all.flatten())      
        distance_sorted = distance_sorted[distance_sorted!=0]
        distance_sorted = distance_sorted[distance_sorted<CN_range]
        
        diff = np.diff(distance_sorted)
        indices = np.where(diff > error_bar)[0] + 1
        self.CN_distances[bond_type] = np.split(distance_sorted, indices)
        self.CN[bond_type] = {np.average(group): group.shape[0]/self.elements_num[center_atom] for group in self.CN_distances[bond_type]}
        if printit is True:
            pprint(self.CN)
        return self.CN
    
    def plot_hist(self, binsize=0.2):
        if not hasattr(self, 'pairs_types'):
            self.get_pairs()
            
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

class ButtonOutputManager:
    def __init__(self):
        self.fc = FileChooser()
        self.fc.register_callback(self.on_file_selected)
        self.output_file = widgets.Output()
        self.fc_box = widgets.HBox([self.fc, self.output_file])
        self.cluster_lst = {}
        
        self._button_layout = widgets.Layout(width='210px', height='25px')
        self._button_layout_narrow = widgets.Layout(width='170px', height='25px')
        self._box_layout = widgets.Layout(display='fixed', 
                                          align_items='center', 
                                          justify_content='flex-start',
                                          )


    def init_menu(self):
        self.output_show = widgets.Output()
        self.output_hist = widgets.Output()
        self.output_getCN = widgets.Output()
        
        self.button_show = widgets.Button(description="Show", layout=self._button_layout_narrow)
        self.button_hist = widgets.Button(description="Plot histogram", layout=self._button_layout_narrow)
        self.button_getCN = widgets.Button(description="Calculate CN", layout=self._button_layout)
        
        self.add_menu_CN()
        self.add_menu_resize()
        self.add_menu_show()

        
        self.menu_box = widgets.HBox([self.box_menu_show,
                                      self.box_menu_CN,
                                      self.box_menu_resize,
                                     ])

        self.output_box = widgets.VBox([self.output_show, 
                                        self.output_hist, 
                                        self.output_getCN])
        
    def update_menu(self):
        options_atoms = [(element, i) for i, element in enumerate(self.cluster.elements_num.keys())]
        self.dropdown_CN_center_atom.options = options_atoms
        self.dropdown_CN_neighbor_atom.options = options_atoms
        self.dropdown_CN_center_atom.value = 0
        self.dropdown_CN_neighbor_atom.value = 0
        self.content_text_shrink_radius.value = self.cluster.get_cluster_size()
        
        self.button_undo_changes.on_click(self.on_button_undo_changes_clicked)
        self.button_show.on_click(self.on_button_show_clicked)
        self.button_hist.on_click(self.on_button_hist_clicked)
        self.button_getCN.on_click(self.on_button_getCN_clicked)
        self.button_shrink.on_click(self.on_button_shrink_clicked)
        self.button_save.on_click(self.on_button_save_clicked)
        self.button_expand.on_click(self.on_button_expand_clicked)
        self.button_remove_atoms.on_click(self.on_button_remove_atoms_clicked)
        self.button_clear.on_click(self.on_button_clear_outpout)

    def add_menu_show(self):
        self.dropdown_view_checkbox_label = widgets.Checkbox(value=True,  description='Labels', disabled=False, layout=widgets.Layout(width='150px'))
        # options_atoms = [(element, i) for i, element in enumerate(self.cluster.elements_num.keys())]
        # self.dropdown_view_highlight_atom1 = widgets.Dropdown(options=options_atoms, value=0, description='Red', layout=widgets.Layout(width='150px'))
        # self.dropdown_view_highlight_atom2 = widgets.Dropdown(options=options_atoms, value=0, description='Blue', layout=widgets.Layout(width='150px'))
        self.content_text_save = widgets.Text(value="Filename", placeholder='File name', description='Save name', disabled=False, layout=widgets.Layout(width='160px'))        
        self.button_save = widgets.Button(description="Save to xyz", layout=self._button_layout_narrow)
        self.button_clear = widgets.Button(description="Clear output", layout=self._button_layout_narrow)

        self.box_menu_show = widgets.VBox([self.button_show,
                                           self.dropdown_view_checkbox_label,
                                           # self.dropdown_view_highlight_atom1,
                                           # self.dropdown_view_highlight_atom2,
                                           self.button_hist, 
                                           self.button_save,
                                           self.content_text_save,
                                           self.button_clear
                                           ])
        
    def add_menu_resize(self):
        self.button_shrink = widgets.Button(description="Shrink cluster", layout=self._button_layout)
        self.content_text_shrink_radius = widgets.FloatText(value=0, placeholder='Radius [A]', description='Radius [A]', disabled=False, step=1, layout=widgets.Layout(width='160px'))
        self.button_expand = widgets.Button(description="Expand CIF", layout=self._button_layout)
        self.content_text_expand_factors_x = widgets.IntText(value=2, placeholder='x', description='x factor', disabled=False, step=1, layout=widgets.Layout(width='160px', align_items='stretch'))
        self.content_text_expand_factors_y = widgets.IntText(value=2, placeholder='y', description='y factor', disabled=False, step=1, layout=widgets.Layout(width='160px'))
        self.content_text_expand_factors_z = widgets.IntText(value=2, placeholder='z', description='z factor', disabled=False, step=1, layout=widgets.Layout(width='160px'))
        # self.content_text_expand_cif = widgets.FloatText(value=0, placeholder='Radius [A]', description='Radius [A]', disabled=False, step=1, layout=widgets.Layout(width='160px'))
        self.button_undo_changes = widgets.Button(description="Undo changes", layout=self._button_layout)
        self.button_remove_atoms = widgets.Button(description="Remove atoms", layout=self._button_layout)
        self.content_text_remove_atoms = widgets.Text(value='', placeholder='e.g. 1, 2, 3', description='atom indices', disabled=False, layout=widgets.Layout(width='160px', align_items='stretch'))
        
        self.box_menu_resize = widgets.VBox([self.button_shrink,
                                             self.content_text_shrink_radius,
                                             self.button_expand,
                                             widgets.VBox([self.content_text_expand_factors_x, 
                                                           self.content_text_expand_factors_y, 
                                                           self.content_text_expand_factors_z], layout=widgets.Layout(width='100%', ajustify_content='flex-start')),
                                             self.button_remove_atoms,
                                             self.content_text_remove_atoms,
                                             self.button_undo_changes,
                                            ])
        
        
    def add_menu_CN(self):
        options = []
        self.dropdown_CN_center_atom = widgets.Dropdown(options=options, value=None, description='Center atom', layout=widgets.Layout(width='180px'))
        self.dropdown_CN_neighbor_atom = widgets.Dropdown(options=options, value=None, description='Neighbor atom', layout=widgets.Layout(width='180px'))
        self.content_text_CN_accuracy = widgets.FloatText(value=0.01, placeholder='Accuracy [A]', description='Accuracy [A]', disabled=False, step=0.01, layout=widgets.Layout(width='180px'))
        self.content_text_CN_range = widgets.FloatText(value=6, placeholder='Range [A]', description='Range [A]', disabled=False, step=1, layout=widgets.Layout(width='180px'))
        
        self.box_menu_CN = widgets.VBox([self.button_getCN, 
                                         self.dropdown_CN_center_atom,
                                         self.dropdown_CN_neighbor_atom,
                                         self.content_text_CN_accuracy,
                                         self.content_text_CN_range
                                        ], layout=self._box_layout)
        
    def init_display(self):
        display(self.fc_box)
        self.init_menu()
        display(self.menu_box)
        display(self.output_box)

    def on_button_remove_atoms_clicked(self, b):
        remove_indices_str = self.content_text_remove_atoms.value
        try:
            remove_indices_lst = [int(i.strip()) for i in remove_indices_str.split(',')]
            # print("Processed list of integers:", remove_indices_lst)
            cluster_atoms_removed = self.cluster.remove_atoms(remove_indices_lst)
            self.cluster.load_xyz(from_file=False, atom_object=cluster_atoms_removed)
            self.on_button_show_clicked(b=None)
            self.content_text_shrink_radius.value = self.cluster.get_cluster_size()
        except:
            with self.output_show:
                clear_output()
                print("Invalid input. Please enter a comma-separated list of numbers.")

    def on_button_show_clicked(self, b):
        with self.output_show:
            clear_output()
            self.cluster.view_xyz(#highlight_atom1=self.dropdown_view_highlight_atom1.label, 
                                  #highlight_atom2=self.dropdown_view_highlight_atom1.label, 
                                  label=self.dropdown_view_checkbox_label.value)
          
    def on_button_hist_clicked(self,b):
        with self.output_hist:
            clear_output()           
            if not hasattr(self, 'pairs_types'):
                self.cluster.get_pairs()
            
            plt.figure(figsize=(8,3))
            for key_i in self.cluster.pairs_group.keys():
                plt.hist(self.cluster.pairs_group[key_i]['distance'], bins=80, alpha=0.3, edgecolor='white', label=key_i)
            plt.xlabel("Distance [Å]")
            plt.ylabel("Number of pairs")
            plt.legend()
            plt.tight_layout()
            plt.show()
            
    def on_button_getCN_clicked(self, b):
        with self.output_getCN:
            clear_output()
            self.cluster.get_CN(center_atom=self.dropdown_CN_center_atom.label, 
                                CN_atom=self.dropdown_CN_neighbor_atom.label, 
                                error_bar=self.content_text_CN_accuracy.value, 
                                CN_range=self.content_text_CN_range.value, 
                                printit=True)
            
    def on_file_selected(self, chooser):
        with self.output_file:
            clear_output()
            print(f"File selected: {chooser.selected}")
            self.cluster = ClusterNeighbor()
            self.cluster.load_xyz(path=self.fc.value)
            self.update_menu()       

    def clear_output_windows(self):
        with self.output_show:
            clear_output()
        with self.output_hist:
            clear_output()
        with self.output_getCN:
            clear_output()  
    
    def on_button_clear_outpout(self, b):
         self.clear_output_windows()
         
    def on_button_undo_changes_clicked(self, b):
        self.cluster.load_xyz(path=self.fc.value)       
        # self.clear_output_windows()
        self.on_button_show_clicked(b=None)
        self.content_text_shrink_radius.value = self.cluster.get_cluster_size()

    
    def on_button_shrink_clicked(self, b):
        new_radius = self.content_text_shrink_radius.value        
        # self.cluster.load_xyz(path=self.fc.value)
        
        if new_radius < self.cluster.get_cluster_size():
            cluster_shrunk = self.cluster.shrink_cluster_size(new_radius=new_radius)
            self.cluster.load_xyz(from_file=False, atom_object=cluster_shrunk)
            # self.clear_output_windows()
            self.on_button_show_clicked(b=None)
        else:
            with self.output_show:
                clear_output()
                print("Should use a radius smaller than cluster size!")
    
    def on_button_expand_clicked(self, b):
        replication_factors = (self.content_text_expand_factors_x.value, 
                               self.content_text_expand_factors_y.value, 
                               self.content_text_expand_factors_z.value)
        try:
            cluster_expanded = self.cluster.expand_cif(replication_factors=replication_factors)
            self.cluster.load_xyz(from_file=False, atom_object=cluster_expanded)
            self.on_button_show_clicked(b=None)
            self.content_text_shrink_radius.value = self.cluster.get_cluster_size()
        except ValueError:
            with self.output_show:
                clear_output()
                print("it's NOT a CIF, cannot expand!")
        
    def on_button_save_clicked(self, b):
        self.cwd = os.path.dirname(self.fc.value)
        filename = os.path.join(self.cwd, self.content_text_save.value + '.xyz')
        print(f"save to {filename}")
        write(filename, self.cluster.atoms)