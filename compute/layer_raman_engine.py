import json
import math
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
from spglib.spglib import get_symmetry_dataset

from .utils.structures import (
    ase_from_tuple,
    get_xsf_structure,
    tuple_from_ase,
    get_covalent_radii_array,
)
from tools_barebone import get_tools_barebone_version
from .utils.hall import hall_numbers_of_spacegroup
from .utils.layers import find_layers, find_common_transformation
from .utils.matrices import matrix_initialization, replace_symbols_with_values
from .utils.pointgroup import (
    pg_number_from_hm_symbol,
    prepare_pointgroup,
    prepare_spacegroup,
    SYMPREC,
)

# Version of this tool
__version__ = "21.11.0"


def nice_print_rot(value, threshold=1.0e-4):
    """
    Converts a float number to a LaTeX string, possibly converting "common" values (integers, and simple square roots)
    to nicer form.

    :param value: a float value
    :param threshold: a numerical threshold to decide if a number is an integer, a square root, ...
    :return: a (LaTeX) string
    """
    int_value = int(round(value))

    if abs(int_value - value) < threshold:
        return f"{int_value:d}"
    if abs(value - 0.5) < threshold:
        return r"\frac{1}{2}"
    if abs(value - (-0.5)) < threshold:
        return r"-\frac{1}{2}"
    if abs(value - math.sqrt(2) / 2) < threshold:
        return r"\frac{\sqrt{2}}{2}"
    if abs(value - (-math.sqrt(2) / 2)) < threshold:
        return r"-\frac{\sqrt{2}}{2}"
    if abs(value - math.sqrt(3) / 2) < threshold:
        return r"\frac{\sqrt{3}}{2}"
    if abs(value - (-math.sqrt(3) / 2)) < threshold:
        return r"-\frac{\sqrt{3}}{2}"

    # As a fallback, return the float representation
    return f"{value:10.5f}"


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

    # Get the primitive cell from the ase cell obtained from the user
    # NOTE! Beside getting the primitive cell, this function will also refine its symmetry.
    primitive_tuple = spglib.find_primitive(
        (
            asecell.get_cell(),
            asecell.get_scaled_positions(),
            asecell.get_atomic_numbers(),
        ),
        symprec=SYMPREC,
    )
    # Get now the conventional cell (it re-does a symmetry analysis)
    dataset = spglib.get_symmetry_dataset(primitive_tuple)
    conventional_tuple = (
        dataset["std_lattice"],
        dataset["std_positions"],
        dataset["std_types"],
    )
    conventional_asecell = ase_from_tuple(conventional_tuple)

    bulk_spg = SpacegroupAnalyzer(
        AseAtomsAdaptor().get_structure(conventional_asecell), symprec=SYMPREC
    )
    pg_bulk_number = pg_number_from_hm_symbol(bulk_spg.get_point_group_symbol())
    return_data["pointgroup_bulk"] = prepare_pointgroup(pg_bulk_number)
    return_data["spacegroup_bulk"] = prepare_spacegroup(bulk_spg)

    # NOTE: there are cases in which it might not be detected - we'll deal with how to display those in the UI

    # From now on, I will work with the conventional cell rather than the one specified by the user
    # This is important because we sometimes (in the output) make assumptions that the number of layers found
    # is the number of layers in the conventional cell (e.g. when we say "Multilayer spacegroup
    # for N >= {num_layers_conventional}").
    is_layered, layer_structures, layer_indices, rotated_asecell = find_layers(
        conventional_asecell, factor=skin_factor
    )

    detected_hall_number = None
    if rotated_asecell is not None:
        # Detect Hall setting
        for hall_number in hall_numbers_of_spacegroup[dataset["number"]]:
            hall_dataset = get_symmetry_dataset(
                tuple_from_ase(rotated_asecell), hall_number=hall_number
            )
            # print(hall_number, hall_dataset['transformation_matrix'], hall_dataset['origin_shift'])

            # If it's Identity, we've identified the correct Hall setting (or at least one among
            # the possible origin choices). We stop at the first one that satisfied this.
            if (
                np.sum(
                    (np.eye(3) - np.array(hall_dataset["transformation_matrix"])) ** 2
                )
                < 1.0e-6
            ):
                detected_hall_number = hall_number
                break

    return_data["hall_number"] = detected_hall_number

    # Get the scaled radii for the bonds detection
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

    rot, transl, center, message = find_common_transformation(
        rotated_asecell, layer_indices
    )

    # Bring back atomic positions so that the origin is the center
    # of the coincidence operation (if found) and atomic positions
    # are inside the the unit cell in the layer plane
    if center is not None:
        rotated_asecell.positions -= center
        rotated_asecell.pbc = [True, True, False]
        rotated_asecell.positions = rotated_asecell.get_positions(wrap=True)
        rotated_asecell.pbc = [True, True, True]
        for layer in layer_structures:
            layer.positions -= center
            layer.pbc = [True, True, False]
            layer.positions = layer.get_positions(wrap=True)
            layer.pbc = [True, True, True]

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

    if rot is None:
        rot_latex = None
    else:
        rot_latex = (
            "R = "
            "\\left(\\begin{array}{ccc}%s & %s & %s \\\\ %s & %s & %s \\\\ %s & %s & %s \\end{array}\\right)"
            % (
                nice_print_rot(rot[0, 0]),
                nice_print_rot(rot[0, 1]),
                nice_print_rot(rot[0, 2]),
                nice_print_rot(rot[1, 0]),
                nice_print_rot(rot[1, 1]),
                nice_print_rot(rot[1, 2]),
                nice_print_rot(rot[2, 0]),
                nice_print_rot(rot[2, 1]),
                nice_print_rot(rot[2, 2]),
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
    return_data["num_layers_bulk"] = num_layers_bulk

    return_data["rotated_cell"] = {
        "layer_cell": rotated_asecell.cell.tolist(),
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

    # As explained also in its docstring, `find_common_transformation` will return a rotation matrix that does NOT
    # flip the z axis, if at least one exists. This is the case for category II (layers are polar, no flip operation)
    # and for category I (layer is non-polar, to it is always possible to combine a coincidence operation with flip
    # together with a bulk operation that also flips z - layers remain invariant and we get rot[2, 2] > 0).
    # If it remains negative, it means it's a category III.
    return_data["is_category_III"] = bool(rot[2, 2] < 0)

    ## COMPUTE HERE VARIOUS POINTGROUP/SPACEGROUP INFORMATION FOR BULK AND VARIOUS MLs
    spg_bilayer = get_symmetry_multilayer(rotated_asecell, layer_indices, num_layers=2)
    all_dicts, all_matrices = construct_all_matrices(
        spg_bilayer,
        num_layers_bulk,
        transformation=rot,
        is_category_III=return_data["is_category_III"],
    )
    fc_dict = construct_force_constant_dict(all_dicts, all_matrices)

    pg_bilayer_number = pg_number_from_hm_symbol(spg_bilayer.get_point_group_symbol())

    spg_monolayer = get_symmetry_multilayer(
        rotated_asecell, layer_indices, num_layers=1
    )
    pg_monolayer_number = pg_number_from_hm_symbol(
        spg_monolayer.get_point_group_symbol()
    )

    # Layer unit-cell surface is the determinant of the 2x2 unit cell
    cell2d = rotated_asecell.cell[:2, :2].tolist()
    layer_surface_ang2 = abs(cell2d[0][0] * cell2d[1][1] - cell2d[0][1] * cell2d[1][0])

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
    num_layers_even = num_layers_bulk + num_layers_bulk % 2
    spg_even = get_symmetry_multilayer(
        rotated_asecell, layer_indices, num_layers=num_layers_even,
    )
    pg_even_number = pg_number_from_hm_symbol(spg_even.get_point_group_symbol())

    # Either a finite ML with num_layers_bulk, or num_layers_bulk + 1 (the odd between the two)
    num_layers_odd = num_layers_bulk + (num_layers_bulk + 1) % 2
    # I need at least two layers to define a multilayer with modes!
    # So if I only have 1, I need to consider 3 as the smallest ML with odd n_primitive
    if num_layers_odd == 1:
        num_layers_odd = 3
    spg_odd = get_symmetry_multilayer(
        rotated_asecell, layer_indices, num_layers=num_layers_odd,
    )
    pg_odd_number = pg_number_from_hm_symbol(spg_odd.get_point_group_symbol())

    # This is possibly different than spg_even for category III - for generality, though,
    # I always compute it
    spg_even_plus_two = get_symmetry_multilayer(
        rotated_asecell, layer_indices, num_layers=num_layers_even + 2
    )
    pg_even_plus_two_number = pg_number_from_hm_symbol(
        spg_even_plus_two.get_point_group_symbol()
    )

    layer_mass_amu = sum(atom.mass for atom in rotated_asecell[layer_indices[0]])

    return_data["pointgroup_even"] = prepare_pointgroup(pg_even_number)
    return_data["pointgroup_odd"] = prepare_pointgroup(pg_odd_number)
    return_data["pointgroup_even_plus_two"] = prepare_pointgroup(
        pg_even_plus_two_number
    )
    # This information is needed, for category III, to decide what to show
    return_data["pointgroup_even_is_multiple_four"] = not (num_layers_even % 4)
    return_data["pointgroup_monolayer"] = prepare_pointgroup(pg_monolayer_number)
    return_data["monolayer_has_z_inversion"] = monolayer_has_z_inversion
    return_data["pointgroup_bilayer"] = prepare_pointgroup(pg_bilayer_number)
    return_data["layer_surface_ang2"] = float(layer_surface_ang2)
    return_data["layer_mass_amu"] = float(layer_mass_amu)

    app_data = {
        "structure": structure,
        "pointgroupEven": pg_even_number,
        "pointgroupOdd": pg_odd_number,
        "uniqueAxisTransformationEven": find_unique_axis_transformation(spg_even),
        "uniqueAxisTransformationOdd": find_unique_axis_transformation(spg_odd),
        "layerSurfaceAng2": float(layer_surface_ang2),
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
        sorted(set({k for m_dict in matrix_dicts for k in m_dict.keys()}))
    ):
        fc_dict["variables"].append(
            {
                "name": v,
                "displayName": string.ascii_lowercase[i],
                "value": matrix_initialization(v),
            }
        )
        # If this happens, we'll probably need to slightly extend the logic
        # to e.g. start using aa, ab, ... but this really only happens for many
        # layers!!
        assert i < 26, "Too many variables!"
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
    cry_sys = spg.get_crystal_system()
    if cry_sys in [
        "tetragonal",
        "hexagonal",
        "trigonal",
    ]:
        matrix_dict = {
            "C11": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            "C33": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        }
    elif cry_sys == "orthorhombic":
        matrix_dict = {
            "C11": [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            "C22": [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            "C33": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        }
    elif cry_sys == "monoclinic":
        matrix_dict = {
            "C11": [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            "C22": [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            "C33": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
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
        matrix_dict.update({"C{}{}".format(*(_ + 1 for _ in nonzero)): matrix})
    elif cry_sys == "triclinic":
        matrix_dict = {
            "C11": [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            "C12": [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            "C13": [[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            "C22": [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            "C23": [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
            "C33": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        }
    elif cry_sys == "cubic":
        raise ValueError(
            "Problems with bilayer structure, cubic crystal system detected"
        )
    else:
        raise ValueError("Bilayer crystal system not valid")
    return matrix_dict


def construct_all_matrices(  # pylint: disable=too-many-locals
    spg, num_layers, transformation, is_category_III, unique_matrix_dict=False
):
    """
    Construct the interlayer force constant matrices given the spacegroup object of the bilayer (N=2)
    and the trasformation matrix that brings EACH layer into the next one.

    Note: in reality only the pointgroup is needed, but for simplicity we pass the whole spacegroup object
    
    :param is_category_III: if True, alternates the matrices at every other layer
       by appending a letter 'a' or 'b' to every constant
    :param unique_matrix_dict: if True, remove duplicates. By default we keep it False, so it's
        easier to explain that the superscript is the layer index the matrix refers to.
    """
    matrix_dict = construct_first_matrix(spg)
    # Under the assumption that this transformation does not flip the z-axis, this
    # is also the transformation that brings a bilayer into the next one.
    matrix_lists = []
    matrix_dicts = []
    this_transformation = np.identity(3)
    for layer_idx in range(num_layers):
        # Alternate matrices for category III
        if is_category_III:
            if layer_idx % 2:
                # 1. we need to use a prefix, since later we check the last characters of the variable
                #    to set some initial random values, and a suffix would break the logic
                # 2. It is good that the string "even" is before "odd" alphabetically, so the symbols
                #    a, b, c, ... remain properly sorted.
                prefix = "odd-"
            else:
                prefix = "even-"
        else:
            prefix = ""
        this_matrix_dict = {
            f"{prefix}{key}": value for key, value in matrix_dict.items()
        }
        m_dict, m_list = rotate_and_simplify_matrix(
            this_matrix_dict, this_transformation
        )
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
    matrix_dict: dictionary whose keys are the free parameters of the interlayer
                  force constant matrix, and the value is a matrix that multiplies the
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
