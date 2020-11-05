import json
import time
import string

import ase
import ase.io
import numpy as np
from ase.data import chemical_symbols
import spglib

from collections.abc import Iterable

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from .utils.structures import (
    ase_from_tuple,
    get_xsf_structure,
    tuple_from_ase,
    get_covalent_radii_array,
)
from tools_barebone import get_tools_barebone_version
from .utils.layers import find_layers, find_common_transformation
from .utils.matrices import matrix_initialization, replace_symbols_with_values
from .utils.pointgroup import (
    pg_number_from_hm_symbol,
    prepare_pointgroup,
    prepare_spacegroup,
    SYMPREC,
)

# Version of this tool
__version__ = "20.11.0"


def process_structure_core(
    structure, logger, flask_request, skin_factor
):  # pylint: disable=unused-argument, too-many-locals, too-many-statements
    start_time = time.time()

    # Get information on the crystal structure to be shown later
    inputstructure_cell_vectors = [
        [idx, coords[0], coords[1], coords[2]]
        for idx, coords in enumerate(structure[0], start=1)
    ]
    inputstructure_symbols = [chemical_symbols[num] for num in structure[2]]
    inputstructure_atoms_scaled = [
        [label, coords[0], coords[1], coords[2]]
        for label, coords in zip(inputstructure_symbols, structure[1])
    ]

    inputstructure_positions_cartesian = np.dot(
        np.array(structure[1]), np.array(structure[0]),
    ).tolist()
    inputstructure_atoms_cartesian = [
        [label, coords[0], coords[1], coords[2]]
        for label, coords in zip(
            inputstructure_symbols, inputstructure_positions_cartesian
        )
    ]

    # prepare template dictionary to return later
    return_data = {
        "app_data_json": json.dumps(
            None
        ),  # None by default, if e.g. layers are not found
        "common_layers_search": None,  # None by default
        "layers": [],  # Empty list if no layers found
        "has_common_layers": False,
        "xsfstructure": get_xsf_structure(structure),
        "inputstructure_cell_vectors": inputstructure_cell_vectors,
        "inputstructure_atoms_scaled": inputstructure_atoms_scaled,
        "inputstructure_atoms_cartesian": inputstructure_atoms_cartesian,
        "skin_factor": skin_factor,
        "spglib_version": spglib.__version__,
        "ase_version": ase.__version__,
        "tools_barebone_version": get_tools_barebone_version(),
        "this_tool_version": __version__,
    }

    asecell = ase_from_tuple(structure)
    is_layered, layer_structures, layer_indices, rotated_asecell = find_layers(
        asecell, factor=skin_factor
    )

    scaled_radii_per_site = skin_factor * get_covalent_radii_array(asecell)
    # This is a dict of the form {"Na": 1.3, "C": 1.5}, ..
    scaled_radii_per_kind = {
        atom.symbol: scaled_radius
        for atom, scaled_radius in zip(asecell, scaled_radii_per_site)
    }

    # I now construct the list of *pairwise* threshold distances, to be passed to JSMOL
    # In theory, I could use simply "set bondTolerance 0;" and "{_P}.bondingRadius = 1.4158" (in this example, for
    # the P element). However, it does not seem to be setting the threshold for showing a bond at the sum,
    # but at some different value.
    # Therefore, I instead compute the pairwise threshold distance, say for elements Ga and As, and pass the following
    # JSMOL string (if, say, I don't want bonds for atoms closer than 0.2 ang, and the threshold distance is 2.27 ang):
    # "connect 0.2 2.27 {_Ga} {_As};"
    # It is good to prepend this with "set autobond off;" before loading, or use first a "connect delete;" to remove
    # existing bonds
    jsmol_bond_commands = []
    min_bonding_distance = 0.2  # ang
    for kind1, radius1 in scaled_radii_per_kind.items():
        for kind2, radius2 in scaled_radii_per_kind.items():
            if kind1 > kind2:
                # Just do one of the two pairs
                continue
            jsmol_bond_commands.append(
                f"connect {min_bonding_distance} {radius1+radius2} {{_{kind1}}} {{_{kind2}}}; "
            )

    # Encode as JSON string before sending, so it's safe to inject in the code
    return_data["jsmol_bond_command"] = json.dumps("".join(jsmol_bond_commands))

    if not is_layered:
        # I return here; some sections will not be present in the output so they will not be shown.
        compute_time = time.time() - start_time
        return_data["compute_time"] = compute_time
        logger.debug(json.dumps(return_data, indent=2, sort_keys=True))
        return return_data

    layer_xsfs = [
        get_xsf_structure(tuple_from_ase(layer_structure))
        for layer_structure in layer_structures
    ]

    return_data["layers"] = list(
        zip(
            layer_xsfs,
            # Needed because this might return int64 numpy objects, not JSON-serializable
            [
                [int(index) for index in this_layer_indices]
                for this_layer_indices in layer_indices
            ],
        )
    )

    rot, transl, message = find_common_transformation(rotated_asecell, layer_indices)

    if rot is None:
        rot_latex = None
    else:
        rot_latex = (
            "R = "
            "\\left(\\begin{array}{ccc}%+10.5f & %+10.5f & %+10.5f \\\\ %+10.5f & %+10.5f & %+10.5f \\\\ %+10.5f & %+10.5f & %+10.5f \\end{array}\\right)"
            % (
                rot[0, 0],
                rot[0, 1],
                rot[0, 2],
                rot[1, 0],
                rot[1, 1],
                rot[1, 2],
                rot[2, 0],
                rot[2, 1],
                rot[2, 2],
            )
        )

    if transl is None:
        transl_latex = None
    else:
        transl_latex = (
            "\\tau = \\left(\\begin{array}{c}%+15.10f \\\\%+15.10f \\\\%+15.10f \\end{array}\\right)\\text{\\AA}"
            % (transl[0], transl[1], transl[2])
        )
    return_data["common_layers_search"] = {
        "rot_latex": rot_latex,
        "transl_latex": transl_latex,
        "message": message,
    }

    # This is returned both in the return_data, for the HTML view, and in the app data,
    # to be set as a minimum for the REST API requests
    num_layers_bulk = len(layer_indices)
    return_data["num_layers_lulk"] = num_layers_bulk

    cell2d = rotated_asecell.cell[:2, :2].tolist()
    return_data["rotated_cell"] = {
        "cell2d": cell2d,
        "layer_atoms": [
            list(
                zip(
                    rotated_asecell[this_layer_indices].symbols,
                    rotated_asecell[this_layer_indices].positions.tolist(),
                )
            )
            for this_layer_indices in layer_indices
        ],
    }

    if rot is None:
        # I return here; some sections will not be present in the output so they will not be shown.
        compute_time = time.time() - start_time
        return_data["compute_time"] = compute_time
        logger.debug(json.dumps(return_data, indent=2, sort_keys=True))
        return return_data

    return_data["has_common_layers"] = True

    ## COMPUTE HERE VARIOUS POINTGROUP/SPACEGROUP INFORMATION FOR BULK AND VARIOUS MLs
    spg_bilayer = get_symmetry_multilayer(rotated_asecell, layer_indices, num_layers=2)
    all_dicts, all_matrices = construct_all_matrices(
        spg_bilayer, num_layers_bulk, transformation=rot
    )
    fc_dict = construct_force_constant_dict(all_dicts, all_matrices)

    pg_bilayer_number = pg_number_from_hm_symbol(spg_bilayer.get_point_group_symbol())

    spg_monolayer = get_symmetry_multilayer(
        rotated_asecell, layer_indices, num_layers=1
    )
    pg_monolayer_number = pg_number_from_hm_symbol(
        spg_monolayer.get_point_group_symbol()
    )

    layer_mass_amu = sum(atom.mass for atom in rotated_asecell[layer_indices[0]])

    # This list contains either -1 or 1; if there is at least one -1, it means that the layer is non-polar
    # (e.g. has inversion symetry or a mirror orthogonal to the stacking axis)
    monolayer_has_z_inversion = any(
        abs(value + 1) < 1.0e-6
        for value in (
            operation.rotation_matrix[2, 2]
            for operation in spg_monolayer.get_point_group_operations()
        )
    )

    # Either a finite ML with num_layers_bulk, or num_layers_bulk + 1 (the even between the two)
    spg_even = get_symmetry_multilayer(
        rotated_asecell,
        layer_indices,
        num_layers=num_layers_bulk + num_layers_bulk % 2,
    )

    pg_even_number = pg_number_from_hm_symbol(spg_even.get_point_group_symbol())
    # Either a finite ML with num_layers_bulk, or num_layers_bulk + 1 (the odd between the two)
    spg_odd = get_symmetry_multilayer(
        rotated_asecell,
        layer_indices,
        num_layers=num_layers_bulk + (num_layers_bulk + 1) % 2,
    )

    pg_odd_number = pg_number_from_hm_symbol(spg_odd.get_point_group_symbol())

    bulk_spg = SpacegroupAnalyzer(
        AseAtomsAdaptor().get_structure(asecell), symprec=SYMPREC
    )
    pg_bulk_number = pg_number_from_hm_symbol(bulk_spg.get_point_group_symbol())

    return_data["pointgroup_even"] = prepare_pointgroup(pg_even_number)
    return_data["pointgroup_odd"] = prepare_pointgroup(pg_odd_number)
    return_data["pointgroup_monolayer"] = prepare_pointgroup(pg_monolayer_number)
    return_data["monolayer_has_z_inversion"] = monolayer_has_z_inversion
    return_data["pointgroup_bilayer"] = prepare_pointgroup(pg_bilayer_number)
    return_data["pointgroup_bulk"] = prepare_pointgroup(pg_bulk_number)
    return_data["spacegroup_bulk"] = prepare_spacegroup(bulk_spg)
    return_data["layer_mass_amu"] = float(layer_mass_amu)
    print(return_data)

    app_data = {
        "structure": structure,
        "pointgroupEven": pg_even_number,
        "pointgroupOdd": pg_odd_number,
        "uniqueAxisTransformationEven": find_unique_axis_transformation(spg_even),
        "uniqueAxisTransformationOdd": find_unique_axis_transformation(spg_odd),
        "layerMassAmu": float(layer_mass_amu),
        # This will be used to decide the mimimum number of layers to show - we don't want to ge below this,
        # as the symmetries might be more.
        "numLayersBulk": num_layers_bulk,
        "forceConstants": fc_dict,
    }
    # Add the JSON to the return_data
    return_data["app_data_json"] = json.dumps(app_data)

    # Add the total compute time and return th dictionary
    compute_time = time.time() - start_time
    return_data["compute_time"] = compute_time
    logger.debug(json.dumps(return_data, indent=2, sort_keys=True))
    return return_data


def construct_force_constant_dict(  # pylint: disable=too-many-locals,too-many-nested-blocks,too-many-locals
    matrix_dicts, matrix_lists
):
    """
    Construct the dictionary with force constant information
    """
    fc_dict = {"matrices": matrix_lists}
    fc_dict.update({"variables": []})
    var_mapping = {}
    for i, v in enumerate(
        sorted({k for m_dict in matrix_dicts for k in m_dict.keys()})
    ):
        fc_dict["variables"].append(
            {
                "name": v,
                "displayName": string.ascii_lowercase[i],
                "value": matrix_initialization(v),
            }
        )
        var_mapping.update({v: string.ascii_lowercase[i]})
    description_list = []
    for ifc, matrix_dict in enumerate(matrix_dicts):
        description = ""
        matrix_dict, matrix = rotate_and_simplify_matrix(matrix_dict, np.identity(3))
        matrix = replace_symbols_with_values(matrix, var_mapping)
        description += "K^{}".format(ifc + 1) + "_{\\alpha\\beta} = "
        description += "\\left(\\begin{array}{ccc}"
        for row in matrix:
            row_text = ""
            for entry_idx, entry in enumerate(row):
                last = entry_idx == len(row) - 1
                if isinstance(entry, Iterable):
                    for variable, factor in entry:
                        if abs(factor - 1.0) < 1e-4:
                            row_text += variable
                        elif abs(factor + 1.0) < 1e-4:
                            row_text += " - " + variable
                        else:
                            row_text += " {:5.3f} ".format(factor) + variable
                else:
                    if abs(entry) < 1e-4:
                        row_text += " 0 "
                    else:
                        row_text += " {:5.3f} ".format(entry)
                if not last:
                    row_text += " & "
            row_text += "\\\\ "
            description += row_text
        description += "\\end{array}\\right)"
        description_list.append(description)
    fc_dict.update({"description": ",\\quad ".join(description_list)})
    return fc_dict


def get_symmetry_multilayer(asecell, layer_indices, num_layers, symprec=SYMPREC):
    """
    Return the SpacegroupAnalyzer object for a finite multilayer with N layers of the given asecell

    IMPORTANT: the layers must be already ordered according to their projection along the stacking direction
    """
    num_layers_bulk = len(layer_indices)
    # define the N-multilayer by taking the first N layers
    # (layers are already ordered according to their projection
    #  along the stacking direction)
    multilayer = asecell[layer_indices[0]]
    for layer_idx in range(1, num_layers):
        z_shift = layer_idx // num_layers_bulk
        new_layer = asecell[layer_indices[layer_idx % num_layers_bulk]]
        new_layer.translate(z_shift * asecell.cell[2])
        multilayer += new_layer

    # put the third lattice of the multiayer orthogonal to the layers
    # and with a large magnitude
    multilayer.cell[2] = [
        0,
        0,
        3.7
        * np.max(
            [
                np.linalg.norm(asecell.cell[0]),
                np.linalg.norm(asecell.cell[1]),
                np.linalg.norm(asecell.cell[2]),
            ]
        ),
    ]

    # transform the structure to a pymatgen structure
    struct = AseAtomsAdaptor().get_structure(multilayer)
    # Find the spacegroup of the multilayer
    spg = SpacegroupAnalyzer(struct, symprec=symprec)

    return spg


def find_unique_axis_transformation(spg):
    """
    Obtain the transformation matrix that brings a unique axis
    along the z direction. Relevant only for monoclinic systems or 
    orthorhombic with pointgroup mm2 (in which case the rotation 
    axis is takes as unique axis)
    For other systems it is either irrelevant (triclinic, other orthorhombic) 
    or already along z by construction (tetragonal, hexagonal, trigonal)
    """
    cry_sys = spg.get_crystal_system()
    if cry_sys == "monoclinic":
        unique_direction = find_unique_axis_monoclinic(spg)
    elif cry_sys == "orthorhombic" and spg.get_point_group_symbol() == "2/m":
        # the unique direction is associated with the two-fold rotation
        for symop in spg.get_point_group_operations(cartesian=True):
            # so we need to discard roto-reflections and the identity
            if np.linalg.det(symop) > 0 and np.sum(np.diag(symop)) < 2.5:
                unique_direction = np.argwhere(np.diag(symop) > 0.0)[0][0]
    else:
        unique_direction = 2

    return rotate_unique_axis(unique_direction).tolist()


def rotate_unique_axis(unique_dir):
    """
    Return the transformation that brings the direction in input
    into the z-direction
    """
    transformation = np.identity(3)
    # if the unique_dir is not already 2 we need to rotate so
    # that the unique_dir becomes 2
    if unique_dir != 2:
        transformation[2, 2] = 0.0
        transformation[unique_dir, unique_dir] = 0.0
        transformation[2, unique_dir] = 1.0
        transformation[unique_dir, 2] = -1.0

    return transformation


def find_unique_axis_monoclinic(spg):
    """
    Find the unique axis for a monoclinic system
    given the pymatgen spacegroup object spg
    """
    # initialize a random tensor
    orig_tensor = np.reshape([np.random.rand() for i in range(9)], (3, 3))
    # impose that the tensor is symmetric
    orig_tensor += orig_tensor.T
    # initialize to zero the symmetrized tensor
    tensor = np.zeros((3, 3))
    # symmetrize tensor by applying all symmetry operations
    for symop in spg.get_point_group_operations(cartesian=True):
        tensor += symop.transform_tensor(orig_tensor)
    # put to zero the diagonal components
    tensor = tensor - np.diag(np.diag(tensor))
    # find non-zero off-diagonal component of the invariant tensor
    nonzero = np.argwhere(abs(tensor) > 1e-5)
    if len(nonzero) != 2:
        raise ValueError(
            "Problems identifying the unique axis in monoclinic system: "
            "more than 2 non-zero off-diagonal entries in the invariant tensor"
        )
    if (nonzero[0] != nonzero[1][::-1]).any():
        raise ValueError(
            "Problems identifying the unique axis in monoclinic system: "
            "invariant tensor not symmetric."
        )
    # The unique axis direction is the one not involved in the off-diagonal matrix element
    direction = list(range(3))
    for i in nonzero[0]:
        direction.remove(i)
    return direction[0]


def construct_first_matrix(spg):
    """
    Construct the interlayer force constant matrix between the first and the second layer
    based on the symmetry of this bilayer. 

    Note: in reality only the pointgroup is needed, but for simplicity we pass the whole spacegroup object
    """
    # TODO: for now it works only when the transformation matrix has no inversion along z
    # (category III, where it is not possible to find one operation without flip along z)
    cry_sys = spg.get_crystal_system()
    if cry_sys in [
        "tetragonal",
        "hexagonal",
        "trigonal",
    ]:
        matrix_dict = {
            "C111": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            "C133": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        }
    elif cry_sys == "orthorhombic":
        matrix_dict = {
            "C111": [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            "C122": [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            "C133": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        }
    elif cry_sys == "monoclinic":
        matrix_dict = {
            "C111": [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            "C122": [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            "C133": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        }
        # Now add the off-diagonal element according to the unique-axis
        # direction in the specific setting
        unique_direction = find_unique_axis_monoclinic(spg)
        nonzero = list(range(3))
        nonzero.remove(unique_direction)
        # Initialize to zero the matrix associated with off-diagonal elements
        matrix = np.zeros((3, 3))
        # and set to one the correct off-diagonal elements
        matrix[[nonzero], [nonzero[::-1]]] = 1.0
        matrix_dict.update({"C1{}{}".format(*(_ + 1 for _ in nonzero)): matrix})
    elif cry_sys == "triclinic":
        matrix_dict = {
            "C111": [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            "C112": [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            "C113": [[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            "C122": [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            "C123": [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
            "C133": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        }
    elif cry_sys == "cubic":
        raise ValueError(
            "Problems with bilayer structure, cubic crystal system detected"
        )
    else:
        raise ValueError("Bilayer crystal system not valid")
    return matrix_dict


def construct_all_matrices(  # pylint: disable=too-many-locals
    spg, num_layers, transformation, unique_matrix_dict=False
):
    """
    Construct the interlayer force constant matrices given the spacegroup object of the bilayer (N=2)
    and the trasformation matrix that brings EACH layer into the next one.

    Note: in reality only the pointgroup is needed, but for simplicity we pass the whole spacegroup object
    """
    matrix_dict = construct_first_matrix(spg)
    # Under the assumption that this transformation does not flip the z-axis, this
    # is also the transformation that brings a bilayer into the next one.
    # TODO: generalize to take into account the case of transformations that flip z
    matrix_lists = []
    matrix_dicts = []
    this_transformation = np.identity(3)
    for _ in range(num_layers):
        m_dict, m_list = rotate_and_simplify_matrix(matrix_dict, this_transformation)
        matrix_lists.append(m_list)
        # In the list of matrix dictionaries add either all of them or check for unicity
        if unique_matrix_dict:
            for orig_dict in matrix_dicts:
                for k, v in orig_dict.items():
                    if np.linalg.norm(np.array(m_dict[k]) - np.array(v)) > 1e-3:
                        matrix_dicts.append(m_dict)
                        break
            if not matrix_dicts:
                matrix_dicts.append(m_dict)
        else:
            matrix_dicts.append(m_dict)
        this_transformation = np.matmul(transformation, this_transformation)
    return matrix_dicts, matrix_lists


def rotate_and_simplify_matrix(matrix_dict, transformation):
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
    new_matrix_dict = {}
    for k, v in matrix_dict.items():
        new_matrix_dict.update(
            {k: np.dot(transformation, np.dot(v, np.linalg.inv(transformation)))}
        )
    # We now convert the dictionary to a list, where each entry is a list of lists of the kind
    # [free parameter, coefficient].
    # If a coefficient is zero the corresponding list is suppressed.
    # If there is only one parameter, we reduce the entry to a simple list of lists
    # (the syntax is: a list of [[p1, a], [p2, b]] denotes the lin. comb. p1 * a + p2 * b)
    # If the list is empty, we replace it with 0.0
    eps = 1e-4  # threshold to consider a coefficient to be 0
    matrix_list = []
    for i in range(3):
        row = []
        for j in range(3):
            entry = []
            for k, v in new_matrix_dict.items():
                if abs(v[i, j]) > eps:
                    entry.append([[k, v[i, j]]])
            if len(entry) == 0:
                row.append(0.0)
            elif len(entry) == 1:
                row.append(entry[0])
            else:
                row.append(entry)
        matrix_list.append(row)
    return new_matrix_dict, matrix_list
