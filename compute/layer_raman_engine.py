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

from .utils.structures import ase_from_tuple, get_xsf_structure, tuple_from_ase
from tools_barebone import get_tools_barebone_version
from .utils.layers import find_layers, find_common_transformation
from .utils.matrices import matrix_initialization, replace_symbols_with_values

# Version of this tool
__version__ = "20.09.0"


def process_structure_core(
    structure, logger, flask_request, skin_factor
):  # pylint: disable=unused-argument, too-many-locals
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
        "app_data_json": None,  # None by default, if e.g. layers are not found
        "common_layers_search": None,  # None by default
        "layers": [],  # Empty list if no layers found
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
            "\\tau = \\left(\\begin{array}{c}%+10.5f \\\\%+10.5f \\\\%+10.5f \\end{array}\\right)\\text{\\AA}"
            % (transl[0], transl[1], transl[2])
        )
    return_data["common_layers_search"] = {
        "rot_latex": rot_latex,
        "transl_latex": transl_latex,
        "message": message,
    }

    all_dicts, all_matrices = construct_all_matrices(
        rotated_asecell, layer_indices, transformation=rot
    )

    fc_dict = construct_force_constant_dict(all_dicts, all_matrices)

    app_data = {
        "structure": structure,
        "symmetryInfo": {},
        "forceConstants": fc_dict,
        #        {
        #            # TO BE FIXED
        #            "description": "K^1_{\\alpha\\beta} = "
        #            "\\left(\\begin{array}{ccc}a & 0 & 0 \\\\ 0 & a & 0 \\\\ 0 & 0 & c \\end{array}\\right), "
        #            "K^2_{\\alpha\\beta} = \\left(\\begin{array}{ccc}a & 0 & 0 \\\\ 0 & a & 0 \\\\ 0 & 0 & c \\end{array}\\right)",
        #            # TO BE FIXED
        #            "variables": [
        #                {"name": "C111", "displayName": "a", "value": 1.0},
        #                {"name": "C133", "displayName": "c", "value": 2.0},
        #            ],
        #            # Check here what should actually be done
        #            "matrices": [all_matrices],
        #        },
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


def construct_all_matrices(  # pylint: disable=too-many-locals
    asecell, layer_indices, transformation, symprec=1e-3
):
    """
    Construct the interlayer force constant matrices from the symmetries 
    of the bilayer and the trasformation matrix that brings EACH layer into
    the next one
    """
    # TODO: for now it works only when the transformation matrix
    # has no inversion along z
    num_layers = len(layer_indices)
    # define the bilayer by taking the first two layers
    # (layers are already ordered according to their projection
    #  along the stacking direction)
    if num_layers > 1:
        bilayer = asecell[layer_indices[0]] + asecell[layer_indices[1]]
    else:
        # If there is only one layer we need to consider the bilayer
        # constructed with its periodic copy
        bilayer = asecell[layer_indices[0]]
        bilayer.translate(asecell.cell[2])
        bilayer += asecell[layer_indices[0]]
    # put the third lattice of the bilayer orthogonal to the layers
    # and with a large magnitude
    bilayer.cell[2] = [
        0,
        0,
        10 * np.max([np.linalg.norm(asecell.cell[0]), np.linalg.norm(asecell.cell[1])]),
    ]
    # transform the structure to a pymatgen structure
    struct = AseAtomsAdaptor().get_structure(bilayer)
    # Find the spacegroup of the bilayer
    spg = SpacegroupAnalyzer(struct, symprec=symprec)
    # print(spg.get_point_group_operations(cartesian=True))
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
        # Initialize to zero the matrix associated with off-diagonal elements
        matrix = np.zeros((3, 3))
        # and set to one the correct off-diagonal elements
        matrix[[nonzero[:, 0]], [nonzero[:, 1]]] = 1.0
        matrix_dict.update({"C1{}{}".format(*(nonzero[0] + 1)): matrix})
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
    # Under the assumption that this transformation does not flip the z-axis, this
    # is also the transformation that brings a bilayer into the next one.
    # TODO: generalize to take into account the case of transformations that flip z
    matrix_lists = []
    matrix_dicts = [matrix_dict]
    this_transformation = np.identity(3)
    for _ in range(num_layers):
        m_dict, m_list = rotate_and_simplify_matrix(matrix_dict, this_transformation)
        matrix_lists.append(m_list)
        for orig_dict in matrix_dicts:
            for k, v in orig_dict.items():
                if np.linalg.norm(np.array(m_dict[k]) - np.array(v)) > 1e-3:
                    matrix_dicts.append(m_dict)
                    break
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
