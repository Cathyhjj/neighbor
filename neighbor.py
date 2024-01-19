# %%
import py3Dmol
import numpy as np
import matplotlib.pyplot as plt

xyz = '''59

O -5.301403 -2.601292 -0.500298
O -3.367214 -1.231646 -1.809521
O -4.379644 0.010097 0.793594
O -2.460089 -3.700719 -0.252618
Pb -3.419867 -1.845311 0.270488
O -3.472519 -2.458976 2.350497
O -1.538331 -1.089330 1.041273
O -2.550761 0.152413 3.644388
O -0.631206 -3.558404 2.598176
Pb -1.590983 -1.702995 3.121282
O -1.643635 -2.316661 5.201291
O 0.290553 -0.947014 3.892068
O -3.346703 0.952956 -3.947334
O -4.359133 2.194699 -1.344219
O -2.424945 3.564345 -2.653442
O -3.437375 4.806088 -0.050327
O -1.517820 1.095271 -1.096539
Pb -2.477597 2.950680 -0.573433
O -2.530250 2.337014 1.506576
O -0.596062 3.706661 0.197353
O -1.608492 4.948403 2.800468
O 0.311064 1.237587 1.754256
Pb -4.068581 1.247684 2.547849
Pb -0.648714 3.092995 2.277362
O -0.701366 2.479330 4.357371
O 1.232822 3.848976 3.048147
O 2.139947 1.379903 4.605050
O -2.139947 -1.379903 -4.605050
O -1.232822 -3.848976 -3.048147
O 0.701366 -2.479330 -4.357371
O -0.311064 -1.237587 -1.754256
O 1.608492 -4.948403 -2.800468
Pb 0.648714 -3.092995 -2.277362
O 0.596062 -3.706661 -0.197353
O 2.530250 -2.337014 -1.506576
O 1.517820 -1.095271 1.096539
O 3.437375 -4.806088 0.050327
Pb -0.942269 -4.795991 0.843921
Pb 2.477597 -2.950680 0.573433
O 2.424945 -3.564345 2.653442
O 4.359133 -2.194699 1.344219
O 3.346703 -0.952956 3.947334
O -0.290553 0.947014 -3.892068
O 1.643635 2.316661 -5.201291
O 0.631206 3.558404 -2.598176
O 2.550761 -0.152413 -3.644388
Pb -1.828884 -0.142316 -2.850795
Pb 1.590983 1.702995 -3.121282
O 1.538331 1.089330 -1.041273
O 3.472519 2.458976 -2.350497
O 2.460089 3.700719 0.252618
O 4.379644 -0.010097 -0.793594
Pb 0.000000 0.000000 0.000000
Pb 3.419867 1.845311 -0.270488
O 3.367214 1.231646 1.809521
O 5.301403 2.601292 0.500298
Pb 1.828884 0.142316 2.850795
Pb 0.942269 4.795991 -0.843921
Pp 4.068581 -1.247684 -2.547849
'''

# %%
compound = xyz.split('\n')[2:-1]
xyz_sorted = [compound[i].split(" ") for i in range(len(compound))]
elements = [xyz_sorted[i][0] for i in range(len(compound))]
coordinates = [np.array([float(xyz_sorted[i][1]),
                         float(xyz_sorted[i][2]), 
                         float(xyz_sorted[i][3])]) for i in range(len(compound))]

# %%
def calculate_distance(coord1, coord2):
    return np.sqrt((coord2[0] - coord1[0])**2 + (coord2[1] - coord1[1])**2 + (coord2[2] - coord1[2])**2)


num_atoms = len(elements)
output = {"distance_all":[],
          "pairs_all":[]}
for i in range(num_atoms):
    atom_i = elements[i]
    for j in range(i + 1, num_atoms):
        atom_j = elements[j]
        distance = calculate_distance(coordinates[i], coordinates[j])
        output['distance_all'].append(distance)
        output['pairs_all'].append(f"{atom_i}-{atom_j}")
        print(f"{atom_i}-{atom_j}: {distance} Å")

# %%
O_O = ['O-O']
Pb_O = ['Pb-O', 'O-Pb', 'Pp-O', 'O-Pp']
Pb_Pb = ['Pb-Pb', 'Pp-Pb']

pairs_sorted = {"O-O":[],
                "Pb-Pb":[],
                "Pb-O":[]}

for index, element in enumerate(output['pairs_all']):
    if element in O_O:
        pairs_sorted["O-O"].append(output['distance_all'][index])
    elif element in Pb_O:
        pairs_sorted["Pb-O"].append(output['distance_all'][index])
    elif element in Pb_Pb:
        pairs_sorted["Pb-Pb"].append(output['distance_all'][index])


# %%
plt.figure()
plt.hist(output['distance_all'], bins=30, color='navy', alpha=0.3,edgecolor='gray')
plt.xlabel("Distance [Å]")
plt.ylabel("Number of pairs")
plt.ylim(0, 180)
plt.figure()


plt.hist(pairs_sorted["O-O"], bins=30, color='red', alpha=0.3, edgecolor='white', label='O-O pairs')
plt.hist(pairs_sorted["Pb-O"], bins=20, color='olive', alpha=0.3,edgecolor='white', label='Pb-O pairs')
plt.hist(pairs_sorted["Pb-Pb"], bins=20, color='purple', alpha=0.5,edgecolor='white',  label='Pb-Pb pairs')
plt.xlabel("Distance [Å]")
plt.ylabel("Number of pairs")
plt.legend()
plt.ylim(0, 180)
plt.show()

# %%
xyzview = py3Dmol.view(width=400,height=400)
xyzview.addModel(xyz,'xyz',)
xyzview.setStyle({'stick':{'radius':.1, 'alpha':0.2, 'color':'gray'}, 
                  'sphere': {'radius':.3}
                  }
                 )

xyzview.addStyle({'atom': 'O'}, 
                 {'sphere': {'color': 'blue', 'radius': 0.3}})  

xyzview.addStyle({'atom': 'Pb'}, 
                 {'sphere': {'color': 'red', 'radius': 0.5}})  

xyzview.addStyle({'atom': 'Pp'}, 
                 {'sphere': {'color': 'yellow', 'radius': 0.5}})  


xyzview.setBackgroundColor('0xeeeeee')
xyzview.zoomTo()
xyzview.show()

# %%
# Create a view object
view = py3Dmol.view(width=800, height=400)
view.addModel(xyz, 'xyz')

# Style the entire molecule
view.setStyle({'stick': {}})

# Change color of a specific atom (e.g., the 5th atom in the file, index 4)
# Note: Py3Dmol's atom indexing starts from 1
view.setStyle({'atomIndex': 1}, {'sphere': {'color': 'red', 'radius': 1}})

# Update the view
view.zoomTo()
view.show()

# %%
import py3Dmol

# Create a view object
view = py3Dmol.view(query='pdb:1ycr')  # Load a molecule, for example, PDB ID 1ycr

# Style the molecule
# view.setStyle({'cartoon': {'color': 'spectrum'}})  # General style for the molecule

# Change color of a specific carbon atom
# Replace 'resi' and 'chain' with the appropriate residue number and chain identifier
# and 'atomIndex' with the index of the carbon atom you want to target
view.addStyle({'atom': 'C', 'resi': 90, 'chain': 'A'}, {'sphere': {'color': 'blue', 'radius': 0.6}})

# Show the molecule
view.show()



# %% [markdown]
# Display local file.

# %%
benz='''
     RDKit          3D

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.9517    0.7811   -0.6622 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2847    1.3329   -0.3121 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2365    0.5518    0.3512 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9517   -0.7811    0.6644 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2847   -1.3329    0.3144 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2365   -0.5518   -0.3489 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END
$$$$'''
view = py3Dmol.view(data=benz,
                    style={'stick':{'colorscheme':'cyanCarbon'}}
                    )
view.show()

# %%
view = py3Dmol.view(query='pdb:1dc9',linked=False,viewergrid=(2,2))
view.setViewStyle({'style':'outline','color':'black','width':0.1})
view.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}},viewer=(0,1))
view.setStyle({'stick':{'colorscheme':'greenCarbon'}},viewer=(1,0))
view.setStyle({'cartoon':{'color':'spectrum'}},viewer=(1,1))
view.removeAllModels(viewer=(0,0))
view.addModel(benz,'sdf',viewer=(0,0))
view.setStyle({'stick':{}},viewer=(0,0))
view.zoomTo(viewer=(0,0))
view.render()


# %%
view = py3Dmol.view(query='pdb:1ycr')
chA = {'chain':'A'}
chB = {'chain':'B'}
view.setStyle(chA,{'cartoon': {'color':'spectrum'}})
view.addSurface(py3Dmol.VDW,{'opacity':0.7,'color':'white'}, chA)
view.setStyle(chB,{'stick':{}})
view.show()

# %%
view = py3Dmol.view(query='pdb:5ire',options={'doAssembly':True})
view.setStyle({'cartoon':{'color':'spectrum'}})
view.show()

# %% [markdown]
# Color by temperature factors

# %%
view = py3Dmol.view(query='pdb:1ycr')
view.setStyle({'cartoon': {'color':'white'}})
view.addSurface(py3Dmol.VDW,{'opacity':0.7,'colorscheme':{'prop':'b','gradient':'sinebow','min':0,'max':70}})

# %%
import requests, base64
r = requests.get('https://mmtf.rcsb.org/v1.0/full/5lgo')
view = py3Dmol.view()
view.addModel(base64.b64encode(r.content).decode(),'mmtf')
view.addUnitCell()
view.zoomTo()

# %% [markdown]
# Specifying individual styles for a viewer grid in the constructor

# %%
view = py3Dmol.view(query='pdb:1dc9',viewergrid=(2,2),style=[[{'stick':{}},{'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}}],
                                                            [{'stick':{'colorscheme':'greenCarbon'}},{'cartoon':{'color':'spectrum'}}]])
view.show()

# %% [markdown]
# It isn't possible to convert Python functions to Javascript functions, but Javascript code can be provided in string form to click/hover callbacks.

# %%
v = py3Dmol.view(query="pdb:1ubq",style={'cartoon':{},'stick':{}})
v.setHoverable({},True,'''function(atom,viewer,event,container) {
                   if(!atom.label) {
                    atom.label = viewer.addLabel(atom.resn+":"+atom.atom,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
                   }}''',
               '''function(atom,viewer) {
                   if(atom.label) {
                    viewer.removeLabel(atom.label);
                    delete atom.label;
                   }
                }''')

# %% [markdown]
# An existing viewer can be modified from different cells with `update`

# %%
benz='''
     RDKit          3D

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.9517    0.7811   -0.6622 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2847    1.3329   -0.3121 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2365    0.5518    0.3512 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9517   -0.7811    0.6644 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2847   -1.3329    0.3144 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2365   -0.5518   -0.3489 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END
$$$$'''
view = py3Dmol.view()
view.addModel(benz,'sdf')
view.setStyle({'stick':{}})
view.zoomTo()
view.show()

# %% [markdown]
# However, **this does not work in colab** because colab sandboxes the JavaScript environments of each cell.

# %%
view.setStyle({'stick':{'color':'blue'}})
view.update()


