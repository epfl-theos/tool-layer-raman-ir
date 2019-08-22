import ase
import numpy as np
import io
import json

class FlaskRedirectException(Exception):
    """
    Class used to return immediately with a flash message and a redirect.
    """
    pass

def parse_structure(filecontent, fileformat, extra_data):
    from .structure_importers import get_structure_tuple, UnknownFormatError

    fileobject = io.StringIO(str(filecontent))
    try:
        structure_tuple = get_structure_tuple(
            fileobject, fileformat, extra_data=extra_data)
    except UnknownFormatError:
        raise FlaskRedirectException("Unknown format '{}'".format(fileformat))
    # Bubble up the exception, will be managed by the level above

    return structure_tuple

def process_structure_core(structure, logger, flask_request, skin_factor):
    """
    Get the data to put in the visualizer jinja2 template
    """

    # This is the data that will be injected in the webpage
    app_data = {
        'structure': None, # Pass info on the structure to display
        'symmetryInfo': {},
        'forceConstants': {
            'description': '\\text{Elastic force constant matrices: }K^1_{\\alpha\\beta} = \\left(\\begin{array}{ccc}a & 0 & 0 \\\\ 0 & a & 0 \\\\ 0 & 0 & c \\end{array}\\right), K^2_{\\alpha\\beta} = \\left(\\begin{array}{ccc}a & 0 & 0 \\\\ 0 & a & 0 \\\\ 0 & 0 & c \\end{array}\\right)', 
            'variables': [
                {
                    'name': 'C111',
                    'displayName': 'a',
                    'value': 1.
                },
                {
                    'name': 'C133',
                    'displayName': 'c',
                    'value': 2.
                },
            ],
            'matrices': [
                [['C111', 0., 0.], [0., 'C111', 0.], [0., 0., 'C133']],
                [[ [['C111', -0.6], ['C133', -0.4]] , 0., 0.], [0., 'C111', 0.], [0., 0., 'C133']]
            ]
        }
    }

    return {
        'test_data': "Some data from the server-side python code: SLIDER: {}; Number of atoms: {}, chemical numbers: {}".format(
            skin_factor, len(structure[1]), structure[2]
        ),
        'app_data_json': json.dumps(app_data)
    }

def _tuple_to_ase(structure_tuple):
    """
    Convert a structure tuple with the format 
    (cell, scaled_positions,element_numbers) to a ASE structure
    
    Currently, if the element numbers do not correspond to a known chemical species,
    an exception is raised. 

    :param structure_tuple: the structure in format (cell, scaled_positions,element_numbers)
    """
    cell, rel_pos, numbers = structure_tuple

    if any([i not in ase.atom.atomic_numbers.values() for i in numbers]):
        raise ValueError
    
    return ase.Atoms(cell=cell,scaled_positions=rel_pos,pbc=True,numbers=numbers)

def get_covalent_radii_array(asecell):
   """ 
   Return a list of the covalent radii for 
   the atoms in the ASE structure using the Cordero values

   :params asecell: the ASE structure
   """
   map_atomic_number_covalent_cordero ={
   1:0.31,  
   2:0.28,
   3:1.28,
   4:0.96,
   5:0.84,
   6:0.76,
   7:0.71,
   8:0.66,
   9:0.57,
   10:0.58,
   11:1.66,
   12:1.41,
   13:1.21,
   14:1.11,
   15:1.07,
   16:1.05,
   17:1.02,
   18:1.06,
   19:2.03,
   20:1.76,
   21:1.7,
   22:1.6,
   23:1.53,
   24:1.39,
   25:1.39,
   26:1.32,
   27:1.26,
   28:1.24,
   29:1.32,
   30:1.22,
   31:1.22,
   32:1.2,
   33:1.19,
   34:1.2,
   35:1.2,
   36:1.16,
   37:2.2,
   38:1.95,
   39:1.9,
   40:1.75,
   41:1.64,
   42:1.54,
   43:1.47,
   44:1.46,
   45:1.42,
   46:1.39,
   47:1.45,
   48:1.44,
   49:1.42,
   50:1.39,
   51:1.39,
   52:1.38,
   53:1.39,
   54:1.4,
   55:2.44,
   56:2.15,
   57:2.07,
   58:2.04,
   59:2.03,
   60:2.01,
   61:1.99,
   62:1.98,
   63:1.98,
   64:1.96,
   65:1.94,
   66:1.92,
   67:1.92,
   68:1.89,
   69:1.9,
   70:1.87,
   71:1.87,
   72:1.75,
   73:1.7,
   74:1.62,
   75:1.51,
   76:1.44,
   77:1.41,
   78:1.36,
   79:1.36,
   80:1.32,
   81:1.45,
   82:1.46,
   83:1.48,
   84:1.4,
   85:1.5,
   86:1.5,
   87:2.6,
   88:2.21,
   89:2.15,
   90:2.06,
   91:2,
   92:1.96,
   93:1.9,
   94:1.87,
   95:1.8,
   96:1.69,
   }
   return np.array([ map_atomic_number_covalent_cordero.get(atom.number) 
            for atom in asecell])

def _shortest_vector_index(array):
    """
    Takes an array of vectors and finds the shortest one.

    :param array: array of vectors
    :return idx: the index of the shortest vector in the array
    """
    idx = np.array([np.linalg.norm(vector) for vector in array]).argmin()
    return idx

def _update_and_rotate_cell(asecell,newcell,layer_indices):
   """
   Update the cell according to the newcell provided,
   and then rotate it so that the first two lattice vectors are in the 
   x-y plane. Atomic positions are refolded moving each layer rigidly.
   """
   from numpy.linalg import norm

   asecell.set_cell(newcell)
   normal_vec = np.cross(newcell[0],newcell[1])
   asecell.rotate(v=normal_vec,a=[0,0,1],center=(0,0,0),rotate_cell=True)
   cell = asecell.cell
   # if the first two lattice vectors have equal magnitude and form
   # a 60deg angle, change the second so that the angle becomes 120
   if (abs(norm(cell[0])-norm(cell[1])) < 1e-6 and
       abs(np.dot(cell[0],cell[1])/np.dot(cell[0],cell[0]) - 0.5) < 1e-3):
      cell[1] -= cell[0]
   asecell.set_cell(cell)
   # finally rotate the first cell vector along x
   asecell.rotate(v=cell[0],a=[1,0,0],center=(0,0,0),rotate_cell=True)
   # Wrap back in the unit cell each layer separately
   for layer in layer_indices:
      # projection of the atomic positions of the layer along the third axis
      proj  = np.dot(asecell.positions[layer],[0,0,1])/asecell.get_volume()
      # move back the vertical position of the layer within the cell
      asecell.positions[layer] -= int(proj.mean()) * asecell.cell[2]
   # fix also the inplane component of the positions   
   asecell.positions[:,:2] = asecell.get_positions(wrap=True)[:,:2]
   # If there is only a single layer, center the atomic positions
   #if len(layer_indices)==1:
   #   asecell.center(axis=2)
   return asecell

def check_neighbors(idx,neighbor_list,asecell,visited,layer):
    """
    Iterative function to get all atoms connected to the idx-th atom

    :param idx: the index of the atom whose neighbors are checked
    :param neighbor_list: the neighbor list object provided by ASE
    :param asecell: the ASE structure
    :param visited: the list of visited sites
    :param layer: a list with atoms belonging to the current layer
    """
    visited.append(idx)
    layer.append(idx)
    indeces, offsets = neighbor_list.get_neighbors(idx)
    for ref, offset in zip(indeces,offsets):
        if ref not in visited:
            if not all(offset ==  np.array([0,0,0])):
                asecell.positions[ref] +=  np.dot(offset,asecell.cell)
                neighbor_list.update(asecell)
            check_neighbors(ref,neighbor_list,asecell,visited,layer)

def _find_layers(asecell,factor=1.1,update_cell=False):
    """
    Obtains all subunits of a given structure by looking
    at the connectivity of the bonds

    """
    from ase.calculators.neighborlist import NeighborList
    from numpy.linalg import norm, matrix_rank

    tol = 1e-6
    nl = NeighborList(factor * get_covalent_radii_array(asecell),
                      bothways=True,self_interaction=False,skin=0.0)
    nl.update(asecell)
    vector1, vector2,vector3 = asecell.cell
    twod_unit = False
    only_twod_units = True
    layer_structures = []
    layer_indices = []
    visited = [] 
    for idx in range(len(asecell)):
       layer = [] 
       if idx not in visited:
          check_neighbors(idx,nl,asecell,visited,layer)
          aselayer = asecell.copy()[layer]
          layer_nl = NeighborList(factor*get_covalent_radii_array(aselayer),
                             bothways=True,self_interaction=False,skin=0.0)
          layer_nl.update(aselayer)
          # We search for the periodic images of the first atom (idx=0)
          # that are connected to at least one atom of the connected layer
          neigh_vec = []
          for idx2 in range(len(aselayer)):
            indices, offsets = layer_nl.get_neighbors(idx2)
            for offset in offsets:
               if not all(offset == [0,0,0]):
                  neigh_vec.append(offset)
          # We define the dimensionality as the rank 
          dim = matrix_rank(neigh_vec)          
          if dim == 2:
             twod_unit = True
             cell = asecell.cell
             vectors = list(np.dot(neigh_vec,cell))
             iv = _shortest_vector_index(vectors)
             vector1 = vectors.pop(iv)
             iv = _shortest_vector_index(vectors)
             vector2 = vectors.pop(iv)
             vector3 = np.cross(vector1,vector2)
             while norm(vector3) < tol:
                iv = _shortest_vector_index(vectors)
                vector2 = vectors.pop(iv)
                vector3 = np.cross(vector1,vector2)
             aselayer = _update_and_rotate_cell(aselayer,[vector1,vector2,vector3],
                                               [range(len(aselayer))])
                  
             disconnected = []
             for i in range(-3,4):
                for j in range(-3,4):
                   for k in range(-3,4):
                      vector = i * cell[0] + j * cell[1] + k * cell[2]
                      if np.dot(vector3,vector) > tol:
                         disconnected.append(vector)
             iv = _shortest_vector_index(disconnected)
             vector3 = disconnected[iv]
          else:
              only_twod_units = False
          layer_structures.append(aselayer)
          layer_indices.append(layer)
      
    if twod_unit and update_cell:
       newcell = [vector1,vector2,vector3]
       #if det(newcell)/det(cell) > 1.0:
          #print "New cell larger than the original cell"
       asecell = _update_and_rotate_cell(asecell,newcell,layer_indices)

    return only_twod_units, asecell[visited], layer_structures

def compare_layers(layers,ltol=0.2,stol=0.3,angle_tol=5.0):
    """
    Compares all layers in the material
    layers:: list of ASE structures corresponding to layers
    ltol:: tolerance on cell length 
    stol:: tolerance on atomic site positions
    angle_tol:: tolerance on cell angles
    """
    from pymatgen.analysis.structure_matcher import StructureMatcher
    from pymatgen.io.ase import AseAtomsAdaptor
    # instance of the adaptor to convert ASE structures to pymatgen format
    adaptor = AseAtomsAdaptor()
    
    # Create an instance of Structure Matcher with the given tolerances
    sm = StructureMatcher(ltol,stol,angle_tol)
    ref_layer = adaptor.get_structure(layers[0])
    all_match = True
    for aselayer in layers[1:]:
        layer = adaptor.get_structure(aselayer)
        if not sm.fit(ref_layer,layer):
            all_match = False
    return all_match