import io
import json

import ase
import ase.io
import numpy as np
from ase.calculators.neighborlist import NeighborList

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.operations import SymmOp
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from .map_positions import are_same_except_order


from tools_barebone.structure_importers import get_structure_tuple, UnknownFormatError


class FlaskRedirectException(Exception):
    """
    Class used to return immediately with a flash message and a redirect.
    """


def parse_structure(filecontent, fileformat, extra_data):
    fileobject = io.StringIO(str(filecontent))
    try:
        structure_tuple = get_structure_tuple(
            fileobject, fileformat, extra_data=extra_data)
    except UnknownFormatError:
        raise FlaskRedirectException("Unknown format '{}'".format(fileformat))
    # Bubble up the exception, will be managed by the level above

    return structure_tuple


def process_structure_core(structure, logger, flask_request, skin_factor):  # pylint: disable=unused-argument
    """
    Get the data to put in the visualizer jinja2 template
    """

    # This is the data that will be injected in the webpage
    app_data = {
        'structure': None, # Pass info on the structure to display
        'symmetryInfo': {},
        'forceConstants': {
            'description': 
                '\\text{Elastic force constant matrices: }K^1_{\\alpha\\beta} = '
                '\\left(\\begin{array}{ccc}a & 0 & 0 \\\\ 0 & a & 0 \\\\ 0 & 0 & c \\end{array}\\right), '
                'K^2_{\\alpha\\beta} = \\left(\\begin{array}{ccc}a & 0 & 0 \\\\ 0 & a & 0 \\\\ 0 & 0 & c \\end{array}\\right)', 
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
    return int(idx)

def _gauss_reduce(vec1,vec2,tol=1e-6):
    """
    Get the shortest vectors in the lattice generated by
    the vectors vec1 and vec2 by using the Gauss reduction method
    """
    reduced = False
    while not reduced:
        length1 = np.linalg.norm(vec1)
        length2 = np.linalg.norm(vec2)
        # First vector should be the shortest between the two
        if (length1 - length2) > tol:
            vec = vec1.copy()
            length = length1
            vec1 = vec2.copy()
            length1 = 1*length2
            vec2 = vec.copy()
            length2 = 1*length
        vec = vec2 - np.round(np.dot(vec1,vec2)/length1**2) * vec1
        length = np.linalg.norm(vec)
        if length1 - length > tol:
            vec2 = vec1.copy()
            vec1 = vec.copy()
        else:
            vec2 = vec.copy()
            reduced = True
    return vec1, vec2

def _update_and_rotate_cell(asecell,newcell,layer_indices):
    """
    Update the cell according to the newcell provided,
    and then rotate it so that the first two lattice vectors are in the 
    x-y plane. Atomic positions are refolded moving each layer rigidly.
    """
    asecell.set_cell(newcell)
    normal_vec = np.cross(newcell[0],newcell[1])
    asecell.rotate(v=normal_vec,a=[0,0,1],center=(0,0,0),rotate_cell=True)
    # it needs to be done twice because of possible bugs in ASE
    normal_vec = np.cross(asecell.cell[0],asecell.cell[1])
    asecell.rotate(v=normal_vec,a=[0,0,1],center=(0,0,0),rotate_cell=True)
    cell = asecell.cell    
    # if the first two lattice vectors have equal magnitude and form
    # a 60deg angle, change the second so that the angle becomes 120
    if ( (abs(np.linalg.norm(cell[0])-np.linalg.norm(cell[1])) < 1e-6) and
         (abs(np.dot(cell[0],cell[1])/np.dot(cell[0],cell[0]) - 0.5) < 1e-3)):
        cell[1] -= cell[0]
    asecell.set_cell(cell)
    # finally rotate the first cell vector along x
    angle = np.arctan2(cell[0,1],cell[0,0])*180/np.pi
    asecell.rotate(-angle,v=[0,0,1],center=(0,0,0),rotate_cell=True)
    # Wrap back in the unit cell each layer separately
    for layer in layer_indices:
        # projection of the atomic positions of the layer along the third axis
        proj  = np.dot(asecell.positions[layer],[0,0,1])
        if len(layer_indices)==1:
            # If there is only a single layer, center the atomic positions
            asecell.positions[layer] -= [0,0,proj.mean()]
        else:
            # move back the vertical position of the layer within the cell
            asecell.positions[layer] -= np.floor(proj.mean()/asecell.cell[2,2]) * asecell.cell[2]
    # fix also the inplane component of the positions   
    asecell.positions[:,:2] = asecell.get_positions(wrap=True)[:,:2]
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

def _find_layers(asecell,factor=1.1,update_cell=True):  # pylint: disable=too-many-locals
    """
    Obtains all subunits of a given structure by looking
    at the connectivity of the bonds

    """
    tol = 1e-6
    nl = NeighborList(factor * get_covalent_radii_array(asecell),
                      bothways=True,self_interaction=False,skin=0.0)
    nl.update(asecell)
    vector1, vector2,vector3 = asecell.cell
    is_layered = True
    layer_structures = []
    layer_indices = []
    visited = [] 
    for idx in range(len(asecell)):  # pylint: disable=too-many-nested-blocks
        layer = [] 
        if idx not in visited:
            check_neighbors(idx,nl,asecell,visited,layer)
            aselayer = asecell.copy()[layer]
            layer_nl = NeighborList(factor * get_covalent_radii_array(aselayer),
                                bothways=True,self_interaction=False,skin=0.0)
            layer_nl.update(aselayer)
            # We search for the periodic images of the first atom (idx=0)
            # that are connected to at least one atom of the connected layer
            neigh_vec = []
            for idx2 in range(len(aselayer)):
                _, offsets = layer_nl.get_neighbors(idx2)
                for offset in offsets:
                    if not all(offset == [0,0,0]):
                        neigh_vec.append(offset)
            # We define the dimensionality as the rank 
            dim = np.linalg.matrix_rank(neigh_vec)          
            if dim == 2:
                cell = asecell.cell
                vectors = list(np.dot(neigh_vec,cell))
                #print vectors
                iv = _shortest_vector_index(vectors)
                vector1 = vectors.pop(iv)
                iv = _shortest_vector_index(vectors)
                vector2 = vectors.pop(iv)
                vector3 = np.cross(vector1,vector2)
                while np.linalg.norm(vector3) < tol:
                    iv = _shortest_vector_index(vectors)
                    vector2 = vectors.pop(iv)
                    vector3 = np.cross(vector1,vector2)
                vector1, vector2 = _gauss_reduce(vector1,vector2)
                vector3 = np.cross(vector1,vector2)
                aselayer = _update_and_rotate_cell(aselayer,[vector1,vector2,vector3],
                                                [list(range(len(aselayer)))])
                disconnected = []
                for i in range(-3,4):
                    for j in range(-3,4):
                        for k in range(-3,4):
                            vector = i * cell[0] + j * cell[1] + k * cell[2]
                            if np.dot(vector3,vector) > tol:
                                disconnected.append(vector)
                iv = _shortest_vector_index(disconnected)
                vector3 = disconnected[iv]
                layer_structures.append(aselayer)
                layer_indices.append(layer)
            else:
                is_layered = False
    if is_layered and update_cell:
        newcell = [vector1,vector2,vector3]
        #print ("BULK")
        #print asecell.cell    
        if abs(np.linalg.det(newcell)/np.linalg.det(cell)-1.0)>1e-3:
            print ("New cell has a different volume the original cell")
        asecell = _update_and_rotate_cell(asecell,newcell,layer_indices)
        # Re-order layers according to their projection 
        # on the stacking direction
        vert_direction = np.cross(asecell.cell[0],asecell.cell[1])
        vert_direction /= np.linalg.norm(vert_direction)
        stack_proj = [ np.dot(layer.positions,vert_direction).mean() 
            for layer in [ asecell[il] for il in layer_indices ] ]
        stack_order = np.argsort(stack_proj)
        # order layers with increasing coordinate along the stacking direction
        layer_indices = [ layer_indices[il] for il in stack_order]
        #print (asecell.cell)
    return is_layered, asecell, layer_indices

def layers_match(layers,ltol=0.2,stol=0.3,angle_tol=5.0):
    """
    Compares all layers in the material
    layers:: list of ASE structures corresponding to layers
    ltol:: tolerance on cell length 
    stol:: tolerance on atomic site positions
    angle_tol:: tolerance on cell angles
    """
    # instance of the adaptor to convert ASE structures to pymatgen format
    adaptor = AseAtomsAdaptor()
    
    # Create an instance of Structure Matcher with the given tolerances
    sm = StructureMatcher(ltol,stol,angle_tol)
    # Translate the refence layer (first) from ASE to pymatgen format
    ref_layer = adaptor.get_structure(layers[0])
    # Initialize to true the variable saying if all layers are identical
    all_match = True
    for aselayer in layers[1:]:
        # Translate layer from ASE to pymatgen format
        layer = adaptor.get_structure(aselayer)
        # If the layer does not match the refence one set to false all_match
        if not sm.fit(ref_layer,layer):
            all_match = False
        #else:
        #    str1, str2, fu, s1_supercell = sm._preprocess(ref_layer,layer)
        #    print (sm._strict_match(str1,str2,fu,s1_supercell,break_on_match=False))
    return all_match

def find_common_transformation(asecell,layer_indices,ltol=0.05,stol=0.05,angle_tol=2.0,symprec=1e-3):  # pylint: disable=too-many-arguments,too-many-locals
    """
    Given an input structure, in ASE format, and the list with 
    the indices of atoms belonging to each layer, determine
    if there exists a common transformation that brings one
    layer to the next

    asecell:: ASE structure of the bulk, where the first two 
              lattice vectors have been re-oriented to lie 
              in the plane of the layers
    layer_indices:: list of lists containing the indices of the
                    atoms belonging to each layer
    ltol:: tolerance on cell length 
    stol:: tolerance on atomic site positions
    angle_tol:: tolerance on cell angles
    """
    # instance of the adaptor to convert ASE structures to pymatgen format
    adaptor = AseAtomsAdaptor()
    
    # Create an instance of Structure Matcher with the given tolerances
    sm = StructureMatcher(ltol,stol,angle_tol,primitive_cell=False)

    # Number of layers
    num_layers = len(layer_indices)
    # If there is only one layer, the transformation is
    # simply a translation along the third axis
    if num_layers == 1:
        return np.eye(3),asecell.cell[2]
    # First check that all layers are identical
    layers = [asecell[layer] for layer in layer_indices]
    if not layers_match(layers):
        # an exception should be raised?
        print ("WARNING: layers are not identical")
        return None, None
    # Layers are already ordered according to their
    # projection along the stacking direction
    # we start by comparing the first and second layer
    layer0 = layers[0]
    layer1 = layers[1]
    str0 = adaptor.get_structure(layer0)
    str1 = adaptor.get_structure(layer1)
    # This is the transformation that brings layer0 into layer1
    transformation01 = sm.get_transformation(str1,str0)
    # define the common lattice of both layers (which is the one of the bulk)
    common_lattice = str0.lattice
    # we use twice the same "common lattice" because they are identical in our case
    # in general we should have str0.lattice and str1.lattice (not necessarily in this order)
    rot01 = np.linalg.solve(np.dot(transformation01[0],common_lattice.matrix),common_lattice.matrix)
    # translation vector in cartesian coordinates
    tr01 = np.matmul(common_lattice.matrix.T,transformation01[1])
    # this is the "symmetry" operation that brings layer0 into layer1
    op01 = SymmOp.from_rotation_and_translation(rot01,tr01)

    spg = SpacegroupAnalyzer(str0,symprec=symprec)
    print(spg.get_symmetry_operations())
    # check that the same operation brings each layer into the next one
    # we need to check not only the symmetry operation found by pymatgen,
    # but also its combination with any of the point group operations of
    # the monolayer: TODO 
    for il in range(1,num_layers):
        layer0 = layers[il]
        layer1 = layers[(il+1)%num_layers]
        # TODO: before transforming, translate back the two layers
        # by il * cell[2]/num_layers
        # the transformed positions of the atoms in the first layer
        pos0 = op01.operate_multi(layer0.positions)
        # that should be identical to the ones of the second layer
        pos1 = layer1.positions
        # we already know from above that the species in each layer are identical
        # so we take the set of the ones in the first layer
        for an in set(layer0.get_chemical_symbols()):
            # check, species by species, that the atomic positions are identical
            # up to a threshold
            # The second variable would be the mapping, but I'm not using it
            distance, _ =  are_same_except_order(
                np.array([ _[0] for _ in zip(pos0,layer0.get_chemical_symbols()) if _[1] == an ]),
                np.array([ _[0] for _ in zip(pos1,layer1.get_chemical_symbols()) if _[1] == an ]),
                common_lattice.matrix)
            # if the distance for any of the atoms of the given species
            # is larger than the threshold the transformation is not the same
            # between all consecutive layers
            if distance.max() > stol:
                # an exception should be raised?
                print ("WARNING: transformation between consecutive layers not always the same")
                return None, None
    return rot01, tr01

def construct_all_matrices(asecell,layer_indices,transformation,symprec=1e-3):
    """
    Construct the interlayer force constant matrices from the symmetries 
    of the bilayer and the trasformation matrix that brings EACH layer into
    the next one
    """
    #TODO: for now it works only when the transformation matrix 
    # has no inversion along z
    # define the bilayer by taking the first two layers
    # (layers are already ordered according to their projection
    #  along the stacking direction)
    bilayer = asecell[layer_indices[0]] + asecell[layer_indices[1]]
    # put the third lattice of the bilayer orthogonal to the layers
    # and with a large magnitude
    bilayer.cell[2] = 10 * np.cross(asecell.cell[0],asecell.cell[1])
    # transform the structure to a pymatgen structure
    struct = AseAtomsAdaptor().get_structure(bilayer)
    # Find the spacegroup of the bilayer
    spg = SpacegroupAnalyzer(struct,symprec=symprec)
    print(spg.get_point_group_symbol(), spg.get_crystal_system())
    cry_sys = spg.get_crystal_system()
    if cry_sys in ['tetragonal','trigonal',]:
        matrix_dict = {'C111': [[1., 0., 0.], [0., 1., 0.], [0., 0., 0.]],
                       'C133': [[0., 0., 0.], [0., 0., 0.], [0., 0., 1.]]
                    }
    elif cry_sys == 'orthorhombic' :
        matrix_dict = {'C111': [[1., 0., 0.], [0., 0., 0.], [0., 0., 0.]],
                       'C122': [[0., 0., 0.], [0., 1., 0.], [0., 0., 0.]],
                       'C133': [[0., 0., 0.], [0., 0., 0.], [0., 0., 1.]]
                    }
    elif cry_sys == 'monoclinic':
        matrix_dict =  {'C111': [[1., 0., 0.], [0., 0., 0.], [0., 0., 0.]],
                       'C122': [[0., 0., 0.], [0., 1., 0.], [0., 0., 0.]],
                       'C133': [[0., 0., 0.], [0., 0., 0.], [0., 0., 1.]]
                    }
        # TODO: add the off-diagonal element according to the unique-axis
        # direction in the specific setting
    elif cry_sys == 'triclinic':
        matrix_dict = {'C111': [[1., 0., 0.], [0., 0., 0.], [0., 0., 0.]],
                       'C112': [[0., 1., 0.], [1., 0., 0.], [0., 0., 0.]],
                       'C113': [[0., 0., 1.], [0., 0., 0.], [1., 0., 0.]],                       
                       'C122': [[0., 0., 0.], [0., 1., 0.], [0., 0., 0.]],
                       'C123': [[0., 0., 0.], [0., 0., 1.], [0., 1., 0.]],
                       'C133': [[0., 0., 0.], [0., 0., 0.], [0., 0., 1.]]
        }
    elif cry_sys == 'cubic':
        raise ValueError("Problems with bilayer structure, cubic crystal system detected")
    else:
        raise ValueError("Bilayer crystal system not valid")
    # Under the assumption that this transformation does not flip the z-axis, this 
    # is also the transformation that brings a bilayer into the next one.
    # TODO: generalize to take into account the case of transformations that flip z
    matrix_list = rotate_and_simplify_matrix(matrix_dict,np.identity(3))
    print(matrix_list)

def rotate_and_simplify_matrix(matrix_dict,transformation):
    """
    This function 'rotates' a matrix written with the internal notation
    matrix_dict:: dictionary whose keys are the free parameters of the interlayer
                  force constant matrix, and the value is a matrix that moltiplies the
                  parameter to obtain the corresponding contribution
    transformation:: 3x3 array that contains the transformation matrix that needs to be applied to 
                     the force constant matrix
    """
    # The new force constant matrix is simply the one obtained by transforming the corresponding
    #  matrix of each free paramater (thanks to the linearity of the operation)
    new_matrix_dict =  {}
    for k, v in matrix_dict.items():
        new_matrix_dict.update({k: np.dot(transformation,np.dot(v,np.linalg.inv(transformation)))})
    # We now convert the dictionary to a list, where each entry is a list of lists of the kind
    # [free parameter, coefficient]. 
    # If a coefficient is zero the corresponding list is suppressed.
    # If there is only one parameter, we reduce the entry to a simple list
    # If the list is empty, we replace it with 0.0
    eps = 1e-4 # threshold to consider a coefficient to be 0
    matrix_list = []
    for i in range(3):
        row = [] 
        for j in range(3):
            entry = []
            for k, v in new_matrix_dict.items():
                if abs(v[i,j]) > eps:
                    entry.append([k,v[i,j]])
            if len(entry)==0:
                row.append(0.0)
            elif len(entry)==1:
                row.append(entry[0])
            else:
                row.append(entry)            
        matrix_list.append(row)
    return matrix_list

def test_structure(filename,factor=1.1):
    asecell = ase.io.read(filename)
    is_layered, asecell, layer_indices = _find_layers(asecell,factor=factor)
    if not is_layered:
        raise ValueError("The material is not layered")
    rot,transl =  find_common_transformation(asecell,layer_indices) 
    print(rot, transl)
    print(construct_all_matrices(asecell,layer_indices,rot))
