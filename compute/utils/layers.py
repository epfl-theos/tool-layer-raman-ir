import numpy as np
from ase.neighborlist import NeighborList
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.operations import SymmOp
from pymatgen.io.ase import AseAtomsAdaptor

from .structures import get_covalent_radii_array
from .linalg import shortest_vector_index, gauss_reduce
from .map_positions import are_same_except_order
from .pointgroup import SYMPREC


def find_layers(  # pylint: disable=too-many-locals,too-many-statements,too-many-branches
    asecell, factor=1.1
):
    """
    Obtains all subunits of a given structure by looking
    at the connectivity of the bonds.

    :param asecell: the bulk unit cell (in ase.Atoms format)
    :param factor: the skin factor
    :return: a tuple with a boolean indicating if the material is layered, a list of layers in the structure (ase format),
        a list of indices of the atoms in each layer, and a rotated bulk ASE cell (with stacking axis along z).
        MOREOVER, it 1) returns layers ordered by stacking index and 2) makes sure the layer is connected when 
        removing the PBC along the third (stacking) axis.
    """
    tol = 1.0e-6
    nl = NeighborList(
        factor * get_covalent_radii_array(asecell),
        bothways=True,
        self_interaction=False,
        skin=0.0,
    )
    nl.update(asecell)
    vector1, vector2, vector3 = asecell.cell
    is_layered = True
    layer_structures = []
    layer_indices = []
    visited = []
    aselayer = None
    final_layered_structures = None

    # Loop over atoms (idx: atom index)
    for idx in range(len(asecell)):  # pylint: disable=too-many-nested-blocks
        # Will contain the indices of the atoms in the "current" layer
        layer = []
        # Check if I already visited this atom
        if idx not in visited:
            # Update 'layer' and 'visited'
            check_neighbors(idx, nl, asecell, visited, layer)
            aselayer = asecell.copy()[layer]
            layer_nl = NeighborList(
                factor * get_covalent_radii_array(aselayer),
                bothways=True,
                self_interaction=False,
                skin=0.0,
            )
            layer_nl.update(aselayer)
            # We search for the periodic images of the first atom (idx=0)
            # that are connected to at least one atom of the connected layer
            neigh_vec = []
            for idx2 in range(len(aselayer)):
                _, offsets = layer_nl.get_neighbors(idx2)
                for offset in offsets:
                    if not all(offset == [0, 0, 0]):
                        neigh_vec.append(offset)
            # We define the dimensionality as the rank
            dim = np.linalg.matrix_rank(neigh_vec)
            if dim == 2:
                cell = asecell.cell
                vectors = list(np.dot(neigh_vec, cell))
                iv = shortest_vector_index(vectors)
                vector1 = vectors.pop(iv)
                iv = shortest_vector_index(vectors)
                vector2 = vectors.pop(iv)
                vector3 = np.cross(vector1, vector2)
                while np.linalg.norm(vector3) < tol:
                    iv = shortest_vector_index(vectors)
                    vector2 = vectors.pop(iv)
                    vector3 = np.cross(vector1, vector2)
                vector1, vector2 = gauss_reduce(vector1, vector2)
                vector3 = np.cross(vector1, vector2)
                aselayer = _update_and_rotate_cell(
                    aselayer, [vector1, vector2, vector3], [list(range(len(aselayer)))]
                )
                disconnected = []
                for i in range(-3, 4):
                    for j in range(-3, 4):
                        for k in range(-3, 4):
                            vector = i * cell[0] + j * cell[1] + k * cell[2]
                            if np.dot(vector3, vector) > tol:
                                disconnected.append(vector)
                iv = shortest_vector_index(disconnected)
                vector3 = disconnected[iv]
                layer_structures.append(aselayer)
                layer_indices.append(layer)
            else:
                is_layered = False
    if is_layered:
        newcell = [vector1, vector2, vector3]
        if abs(np.linalg.det(newcell) / np.linalg.det(cell) - 1.0) > 1e-3:
            raise ValueError(
                "An error occurred. The new cell after rotation has a different volume than the original cell"
            )
        rotated_asecell = _update_and_rotate_cell(asecell, newcell, layer_indices)
        # Re-order layers according to their projection
        # on the stacking direction
        vert_direction = np.cross(rotated_asecell.cell[0], rotated_asecell.cell[1])
        vert_direction /= np.linalg.norm(vert_direction)
        stack_proj = [
            np.dot(layer.positions, vert_direction).mean()
            for layer in [rotated_asecell[il] for il in layer_indices]
        ]
        stack_order = np.argsort(stack_proj)
        # order layers with increasing coordinate along the stacking direction
        layer_indices = [layer_indices[il] for il in stack_order]

        # Move the atoms along the third lattice vector so that
        # the first layer has zero projection along the vertical direction
        trans_vector = -(
            stack_proj[stack_order[0]]
            / np.dot(vert_direction, rotated_asecell.cell[2])
            * rotated_asecell.cell[2]
        )
        rotated_asecell.translate(trans_vector)

        # I don't return the 'layer_structures' because there the atoms are moved
        # from their positions and the z axis lenght might not be appropriate
        final_layered_structures = [
            rotated_asecell[this_layer_indices] for this_layer_indices in layer_indices
        ]
    else:
        rotated_asecell = None

    if not is_layered:
        aselayer = None
    return is_layered, final_layered_structures, layer_indices, rotated_asecell


def layers_match(layers, ltol=0.2, stol=0.3, angle_tol=5.0):
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
    sm = StructureMatcher(ltol, stol, angle_tol)
    # Translate the refence layer (first) from ASE to pymatgen format
    ref_layer = adaptor.get_structure(layers[0])
    # Initialize to true the variable saying if all layers are identical
    all_match = True
    for aselayer in layers[1:]:
        # Translate layer from ASE to pymatgen format
        layer = adaptor.get_structure(aselayer)
        # If the layer does not match the refence one set to false all_match
        if not sm.fit(ref_layer, layer):
            all_match = False
        # else:
        #    str1, str2, fu, s1_supercell = sm._preprocess(ref_layer,layer)
        #    print (sm._strict_match(str1,str2,fu,s1_supercell,break_on_match=False))
    return all_match


def check_neighbors(idx, neighbor_list, asecell, visited, layer):
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
    for ref, offset in zip(indeces, offsets):
        if ref not in visited:
            if not all(offset == np.array([0, 0, 0])):
                asecell.positions[ref] += np.dot(offset, asecell.cell)
                neighbor_list.update(asecell)
            check_neighbors(ref, neighbor_list, asecell, visited, layer)


def _update_and_rotate_cell(asecell, newcell, layer_indices):
    """
    Update the cell according to the newcell provided,
    and then rotate it so that the first two lattice vectors are in the 
    x-y plane. Atomic positions are refolded moving each layer rigidly.
    """
    asecell.set_cell(newcell)
    normal_vec = np.cross(newcell[0], newcell[1])
    asecell.rotate(v=normal_vec, a=[0, 0, 1], center=(0, 0, 0), rotate_cell=True)
    # it needs to be done twice because of possible bugs in ASE
    normal_vec = np.cross(asecell.cell[0], asecell.cell[1])
    asecell.rotate(v=normal_vec, a=[0, 0, 1], center=(0, 0, 0), rotate_cell=True)
    cell = asecell.cell
    # if the first two lattice vectors have equal magnitude and form
    # a 60deg angle, change the second so that the angle becomes 120
    if (abs(np.linalg.norm(cell[0]) - np.linalg.norm(cell[1])) < 1e-6) and (
        abs(np.dot(cell[0], cell[1]) / np.dot(cell[0], cell[0]) - 0.5) < 1e-3
    ):
        cell[1] -= cell[0]
    asecell.set_cell(cell)
    # finally rotate the first cell vector along x
    angle = np.arctan2(cell[0, 1], cell[0, 0]) * 180 / np.pi
    asecell.rotate(-angle, v=[0, 0, 1], center=(0, 0, 0), rotate_cell=True)
    # Wrap back in the unit cell each layer separately
    for layer in layer_indices:
        # projection of the atomic positions of the layer along the third axis
        proj = np.dot(asecell.positions[layer], [0, 0, 1])
        if len(layer_indices) == 1:
            # If there is only a single layer, center the atomic positions
            asecell.positions[layer] -= (
                proj.mean() / asecell.cell[2, 2] * asecell.cell[2]
            )
        else:
            # move back the vertical position of the layer within the cell
            asecell.positions[layer] -= (
                np.floor(proj.mean() / asecell.cell[2, 2]) * asecell.cell[2]
            )
    # fix also the inplane component of the positions
    asecell.pbc = [True, True, False]
    asecell.positions = asecell.get_positions(wrap=True)
    asecell.pbc = [True, True, True]
    return asecell


def find_common_transformation(
    asecell, layer_indices, ltol=0.05, stol=0.05, angle_tol=2.0
):  # pylint: disable=too-many-arguments,too-many-locals,too-many-branches
    """
    Given an input structure, in ASE format, and the list with 
    the indices of atoms belonging to each layer, determine
    if there exists a common transformation that brings one
    layer to the next.

    :note: it MUST receive the layers already rearranged so that they are ordered
        along the stacking axis, and with all atoms nearby (i.e., when removing the PBC
        along the third axis, I should still see the layer and not it broken in two or more parts).
        This is performed by the function `find_layers`.

    :param asecell: ASE structure of the bulk, where the first two 
              lattice vectors have been re-oriented to lie 
              in the plane of the layers
    :param layer_indices: list of lists containing the indices of the
                    atoms belonging to each layer
    :param ltol: tolerance on cell length 
    :param stol: tolerance on atomic site positions
    :param angle_tol: tolerance on cell angles
    :return: a tuple of length three: either (rot, transl, None) if there is a common transformation,
        or (None, None, message) if a common transformation could not be found.
        IMPORTANT: the code will try to find a transformation matrix so that rot[2,2] > 0, if it exists.
        If it does not find it, and it returns rot[2,2] < 0, it means it's a category III material (see
        definition in the text). Otherwise, it's either category I or II depending on the monolayer symmetry).
    """
    # instance of the adaptor to convert ASE structures to pymatgen format
    adaptor = AseAtomsAdaptor()

    # Create an instance of Structure Matcher with the given tolerances
    sm = StructureMatcher(ltol, stol, angle_tol, primitive_cell=False)

    # Number of layers
    num_layers = len(layer_indices)
    # If there is only one layer, the transformation is
    # simply a translation along the third axis
    if num_layers == 1:
        return np.eye(3), asecell.cell[2], None

    # If we are here, there are at least two layers
    # First check that all layers are identical
    layers = [asecell[layer] for layer in layer_indices]
    if not layers_match(layers, stol=stol):
        return None, None, "Layers are not identical"

    # Layers are already ordered according to their
    # projection along the stacking direction
    # we start by comparing the first and second layer
    layer0 = layers[0].copy()
    layer0.cell[2] = [0, 0, asecell.cell[2, 2]]
    layer1 = layers[1].copy()
    layer1.cell[2] = [0, 0, asecell.cell[2, 2]]
    str0 = adaptor.get_structure(layer0)
    str1 = adaptor.get_structure(layer1)
    # This is the transformation that brings layer0 onto layer1
    # it is a tuple with a supercell matrix, a translation vector, and the indices
    # of how atoms are rearranged
    transformation01 = sm.get_transformation(str1, str0)

    # There are still cases that, even if layers_match was OK, here we get `None` for transformation.
    # Probably the difference with the code above is in the thresholds and/or in the way the pymatgen
    # code internally implements it. Note: probably we could skip the call to `layers_match` above
    # and just keep the logic below!
    # To avoid crashes, we cope with this here
    if transformation01 is None:
        return None, None, "Layers are not identical"

    # define the common lattice of both layers (which is the one of the bulk)
    # we use twice the same "common lattice" because they are identical in our case
    # in general we should have str0.lattice and str1.lattice (not necessarily in this order)
    common_lattice = str0.lattice
    # find the expression of the rotation in the common lattice
    rot01 = np.linalg.solve(
        np.dot(transformation01[0], common_lattice.matrix), common_lattice.matrix
    )
    # translation vector in Cartesian coordinates
    tr01 = np.matmul(common_lattice.matrix.T, transformation01[1])
    # this is the "symmetry" operation that brings layer0 into layer1
    op01 = SymmOp.from_rotation_and_translation(rot01, tr01)

    # While on the x-y plane we are OK if we find a periodic image, on z we want to
    # drop periodic boundary conditions, so finding a periodic image (using the periodicity
    # of the bulk) is not OK.
    # We therefore check if we can identify the same translation along z for all atoms
    # Note that I only take the z coordinate (last column); this is in angstrom
    required_z_shift_per_atom = (
        layers[1].positions
        - op01.operate_multi(layers[0].positions)[transformation01[2]]
    )[:, 2]
    # Note however that if the layers are not properly defined with all atoms closeby,
    # The values of this array might not be all identical.
    # E.g. if atoms of the same layer are not closeby unless one considers PBC.

    z_shift = required_z_shift_per_atom[0]
    for atom_idx in range(1, len(required_z_shift_per_atom)):
        # NOTE: This test is going to work because the function called before this
        # already reorganised layers to have all atoms close to each other
        assert np.abs(z_shift - required_z_shift_per_atom[atom_idx]) < 5.0e-3

    # Get the spacegroup of the monolayer (useful below)
    spg_mono = SpacegroupAnalyzer(adaptor.get_structure(layer0), symprec=SYMPREC)
    # If the transformation involves a flip in the z-direction
    # we check if, by combining it with a symmetry of the monolayer,
    # we get a transformation that does NOT flip z
    # As soon as I find one, I replace tr01, rot01 and op01 (the combination of tr01 and rot01)
    # with those we found, so that the axis is not flipped and transformation01[0][2, 2] > 0
    if transformation01[0][2, 2] < 0:
        for op in spg_mono.get_symmetry_operations(cartesian=True):
            affine_prod = np.dot(op01.affine_matrix, op.affine_matrix)
            if affine_prod[2, 2] > 0:
                tr01 = affine_prod[0:3][:, 3]
                rot01 = affine_prod[0:3][:, 0:3]
                op01 = SymmOp.from_rotation_and_translation(rot01, tr01)
                break

    # Here, there are now two options:
    # - I couldn't find any operation: then I still have transformation01[0][2, 2] < 0
    #   and I will categorize it as category III
    # - I could find it, then transformation01[0][2, 2] > 0 and it's either category I or II
    #   (depending on the monolayer symmetry)
    #   It could also be category III in weird cases like identical layers that invariant
    #   under z -> -z but with two alternating interlayer distances, i.e. "dimerized"

    # I now use op01 to check if, with op01, I can bring each layer onto the next,
    # and layer num_layers onto num_layers+1, i.e. the first one + the third lattice vector
    # If op01 does not work we might need to combine it with symmetry operations of the monolayer
    # before concluding that the system is not a MDO polytype

    print(op01)

    for symop in spg_mono.get_symmetry_operations(cartesian=True):
        # symmetry operations of the monolayer that flip z are possible
        # only in category I, but would result in a coincidence operation
        # that flip z, which is not necessary in this case (category I), as in this case
        # we'll find another one that does the same job without flipping, so we skip these operations.
        # (also because we want to give the guarantee that, if symop.affine_matrix[2, 2] < 0, then we are
        # in category III)
        if symop.affine_matrix[2, 2] < 0:
            continue
        # combine the coincidence operation with the
        affine_prod = np.dot(op01.affine_matrix, symop.affine_matrix)
        this_rot = affine_prod[0:3][:, 0:3]
        this_tr = affine_prod[0:3][:, 3]
        # To be compatible with the bulk periodicity in the third vertical direction
        # the fractional translation should be such that, multiplied by the number of layers,
        # it gives a Bravais lattice vector
        # If this is not the case it means that the fractional translation has a contribution
        # arising from the fact that the operation has to carried out around a point which is not the
        # origin
        # Let's save in vec the remainder of num_layers * this_tr modulo a Bravais lattice vector,
        # then divided by the number of layers
        vec = (
            np.dot(
                (np.dot(num_layers * this_tr + 1e-12, np.linalg.inv(asecell.cell))) % 1,
                asecell.cell,
            )
            / num_layers
        )
        # If we subtract vec from this_tr, then num_layers * this_tr is a Bravais lattice vector
        this_tr -= vec
        # If vec is non zero it means that the coincidence operation needs to be carried out around a point
        # with coordinates this_tr0 which is not the origin. In particular, the relationship between the two is
        # vec = this_tr0 - this_rot * this_tr0
        # that we need to invert using a pseudoinverse to get this_tr0
        if np.linalg.norm(vec) > 1e-4:
            this_tr0 = np.dot(np.linalg.pinv(np.identity(3) - this_rot), vec)
        else:
            this_tr0 = np.array([0, 0, 0])
        # The coincidence operation is then defined in terms of this_tr only
        this_op = SymmOp.from_rotation_and_translation(this_rot, this_tr)

        found_common = True
        # NOTE: here we start from layer ZERO! The reason is that we want to
        # double check that now that we made a combination with a symmetry operation,
        # we still send layer 1 into layer 2.
        # So we want to check that case as well.
        for il in range(0, num_layers):
            # We need to copy the layers as we'll change them in place
            layer0 = layers[il].copy()
            layer1 = layers[(il + 1) % num_layers].copy()
            # translate back the two layers to put at the origin the
            # center around which we want to perform the coincidence operation
            # if layer1 is the layer num_layer + 1 we need to translate it
            # by a full lattice vector
            layer0.translate(-il * this_tr - this_tr0)
            layer1.translate(
                (
                    -il * this_tr
                    - this_tr0
                    + np.floor((il + 1.0) / num_layers) * asecell.cell[2]
                )
            )
            # the transformed positions of the atoms in the first layer
            pos0 = this_op.operate_multi(layer0.positions)
            # that should be identical to the ones of the second layer
            pos1 = layer1.positions
            # we already know from above that the species in each layer are identical
            # so we take the set of the ones in the first layer
            for an in set(layer0.get_chemical_symbols()):
                # check, species by species, that the atomic positions are identical
                # up to a threshold
                # The second variable would be the mapping, but I'm not using it
                distance, _ = are_same_except_order(
                    np.array(
                        [
                            _[0]
                            for _ in zip(pos0, layer0.get_chemical_symbols())
                            if _[1] == an
                        ]
                    ),
                    np.array(
                        [
                            _[0]
                            for _ in zip(pos1, layer1.get_chemical_symbols())
                            if _[1] == an
                        ]
                    ),
                    common_lattice.matrix,
                )
                # if the distance for any of the atoms of the given species
                # is larger than the threshold the transformation is not the same
                # between all consecutive layers
                found_common = found_common and distance.max() < stol
                # if this transformation does not work it is useless con continue with
                # other species of this layer
                if not found_common:
                    break
            # if this transformation does not work it is useless con continue with
            # subsequent layers
            if not found_common:
                break
        # if the transformation works, no need to test additional transformations
        if found_common:
            break
    if not found_common:
        return (
            None,
            None,
            "The transformation between consecutive layers is not always the same",
        )
    return this_rot, this_tr, None
