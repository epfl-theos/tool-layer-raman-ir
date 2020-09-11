import json
import logging
import os

import flask
import numpy as np

from collections.abc import Iterable

from .layer_raman_engine import (
    parse_structure,
    process_structure_core,
    FlaskRedirectException,
)

verbose = True

logger = logging.getLogger("layer-raman-tool-app")
blueprint = flask.Blueprint("compute", __name__, url_prefix="/compute")


@blueprint.route("/process_structure/", methods=["GET", "POST"])
def process_structure():
    if flask.request.method == "POST":
        # check if the post request has the file part
        if "structurefile" not in flask.request.files:
            return flask.redirect(flask.url_for("input_structure"))
        structurefile = flask.request.files["structurefile"]
        fileformat = flask.request.form.get("fileformat", "unknown")
        filecontent = structurefile.read().decode("utf-8")
        skin_factor = flask.request.form.get("skin-factor", "")
        try:
            skin_factor = float(skin_factor)
        except ValueError:
            flask.flash(
                "Invalid value for skin factor, must be float ({})".format(skin_factor)
            )
            return flask.redirect(flask.url_for("input_data"))

        try:
            structure = parse_structure(
                filecontent=filecontent,
                fileformat=fileformat,
                extra_data=dict(flask.request.form),
            )
        except Exception as exc:
            flask.flash(
                "Unable to parse the structure, sorry... ({}, {})".format(
                    str(type(exc)), str(exc)
                )
            )
            return flask.redirect(flask.url_for("input_data"))

        try:
            data_for_template = process_structure_core(
                structure=structure,
                logger=logger,
                flask_request=flask.request,
                skin_factor=skin_factor,
            )
            return flask.render_template(
                "user_templates/visualizer.html", **data_for_template
            )
        except FlaskRedirectException as e:
            flask.flash(str(e))
            return flask.redirect(flask.url_for("input_data"))
        except Exception as exc:
            flask.flash(
                "Unable to process the structure, sorry... ({}, {})".format(
                    str(type(exc)), str(exc)
                )
            )
            return flask.redirect(flask.url_for("input_data"))
    else:  # GET Request
        return flask.redirect(flask.url_for("input_data"))


VALID_EXAMPLES = {
    # MoS2 bulk from Materials Cloud:
    # https://www.materialscloud.org/explore/2dstructures/details/6e58409f-4ab2-4883-9686-87d4d89c0bf9
    # (Originally from COD, 9007660, P6_3/mmc)
    "MoS2": ("6e58409f-4ab2-4883-9686-87d4d89c0bf9_structure.xsf", 1.1),
}


@blueprint.route("/process_example_structure/", methods=["GET", "POST"])
def process_example_structure():
    """
    Process an example structure (example name from POST request)
    """
    if flask.request.method == "POST":
        examplestructure = flask.request.form.get("examplestructure", "<none>")
        fileformat = "xsf-ase"

        try:
            filename, skin_factor = VALID_EXAMPLES[examplestructure]
        except KeyError:
            flask.flash("Invalid example structure '{}'".format(examplestructure))
            return flask.redirect(flask.url_for("input_data"))

        # I expect that the valid_examples dictionary already filters only
        # existing files, so I don't try/except here
        with open(
            os.path.join(os.path.dirname(__file__), "xsf-examples", filename,)
        ) as structurefile:
            filecontent = structurefile.read()

        try:
            structure = parse_structure(filecontent=filecontent, fileformat=fileformat,)
        except Exception as exc:
            flask.flash(
                "Unable to parse the example structure, sorry... ({}, {})".format(
                    str(type(exc)), str(exc)
                )
            )
            return flask.redirect(flask.url_for("input_data"))

        try:
            data_for_template = process_structure_core(
                structure=structure,
                logger=logger,
                flask_request=flask.request,
                skin_factor=skin_factor,
            )
            return flask.render_template(
                "user_templates/visualizer.html", **data_for_template
            )
        except FlaskRedirectException as e:
            flask.flash(str(e))
            return flask.redirect(flask.url_for("input_data"))
        except Exception as exc:
            flask.flash(
                "Unable to process the structure, sorry... ({}, {})".format(
                    str(type(exc)), str(exc)
                )
            )
            return flask.redirect(flask.url_for("input_data"))
    else:  # GET Request
        return flask.redirect(flask.url_for("input_data"))


def make_response(message, error_code):
    if verbose:
        print(message)
    return flask.make_response(message, error_code)


def get_block_coordinates(i, j, block_size=3):
    """Return numpy slice coordinates for block at (i, j) *block* coordinates.
    
    :param i: row block coordinate (for a block of size `block_size x block_size`)
    :param j: column block coordinate (for a block of size `block_size x block_size`)
    :param block_size: the block size
    """
    return (
        slice(block_size * i, block_size * (i + 1)),
        slice(block_size * j, block_size * (j + 1)),
    )


@blueprint.route("/api/modes/", methods=["POST"])
def get_modes():
    if flask.request.method != "POST":
        return make_response("Only POST method allowed", 405)

    # If not present, assume JSON, even if by default it is usually application/x-www-form-urlencoded
    content_type = flask.request.headers.get("Content-Type", "application/json")
    if content_type.partition(";")[0].strip() != "application/json":
        print(content_type.partition(";")[0].strip())
        return make_response(
            "Invalid request, Content-Type must be 'application/json' instead of '{}'".format(
                content_type
            ),
            400,
        )

    try:
        data = json.loads(flask.request.data.decode("utf-8"))
    except ValueError:
        print(flask.request.data.decode("utf-8"))
        return make_response("Invalid request, not a valid JSON", 400)

    try:
        max_layers = int(data["maxLayers"])
        if max_layers < 2 or max_layers > 40:
            raise KeyError
    except (KeyError, ValueError):
        return make_response(
            "Invalid request, maxLayers value not passed or not in valid range", 400
        )

    try:
        force_constant_params = data["forceConstantParams"]
    except KeyError:
        return make_response("Invalid request, missing forceConstantParams", 400)

    try:
        symbolic_matrices = data["matrices"]
    except KeyError:
        return make_response("Invalid request, missing matrices", 400)

    ## an example of validation - to decide if we want to do it!
    # if 'C111' not in force_constant_params or force_constant_params['C111'] < 0:
    #    return make_response("missing or invalid C111", 400)

    numeric_matrices = replace_symbols_with_values(
        list_of_lists=symbolic_matrices, replacements=force_constant_params
    )
    ## This still contains information on linear combinations of coefficients.
    ## Now, we assume it's a list of 3x3 matrices, and if an element
    ## is not a numeric value but a list of lists, it means a linear combination:
    ## replace with correct combination. E.g.:
    ## [['C111', -1]] = -C111;
    ## [['C111', -0.5], ['C112', -0.5]] = -0.5 * C111 - 0.5 * C112
    numeric_matrices = replace_linear_combinations(numeric_matrices)

    # print(symbolic_matrices)
    # print(force_constant_params)
    # print(numeric_matrices)

    ## LOGIC START ##

    # TODO: check, this is just a first attempt

    # I have to solve the eq. of motion - omega^2 U = K U

    plot_data_x = []
    plot_data_y = []
    # TODO: add layer mass here should be a constant for all layers! and fix units
    for num_layers in range(2, max_layers + 1):
        K_matrix = np.zeros((num_layers * 3, num_layers * 3))

        for block_idx in range(num_layers):
            # Up interaction
            if block_idx < num_layers - 1:  # Not in the last layer
                current_block = numeric_matrices[block_idx % len(numeric_matrices)]
                K_matrix[get_block_coordinates(block_idx, block_idx)] += current_block
                K_matrix[
                    get_block_coordinates(block_idx + 1, block_idx)
                ] -= current_block
            # Down interaction
            if block_idx > 0:  # Not in the first layer
                previous_block = numeric_matrices[
                    (block_idx - 1) % len(numeric_matrices)
                ]
                K_matrix[get_block_coordinates(block_idx, block_idx)] += previous_block
                K_matrix[
                    get_block_coordinates(block_idx - 1, block_idx)
                ] -= previous_block

        eigvals, eigvects = np.linalg.eigh(K_matrix)

        # The first three should be acoustic i.e. almost zero; the rest should be positive
        assert np.sum(np.abs(eigvals[:3])) < 1.0e-10

        plot_data_x += [num_layers] * 3 * (num_layers - 1)
        plot_data_y += eigvals[3:].tolist()

    return_data = {
        "x": plot_data_x,
        "y": plot_data_y,
        # TODO: THIS IS RANDOM!!! TO BE FIXED
        "isBackScattering": [
            idx % 2 < 1 and (idx % 4 < 2) for idx in range(len(plot_data_x))
        ],  # If it's not Raman active, then this is always False
        "isRamanActive": [idx % 4 < 2 for idx in range(len(plot_data_x))],
        "isInfraredActive": [idx % 8 < 4 for idx in range(len(plot_data_x))],
    }

    ## LOGIC END ##

    return flask.jsonify(return_data)


def replace_symbols_with_values(list_of_lists, replacements):
    """
    Take iteratively a list of lists (at any depth level) and replace strings with
    the corresponding values in the ``replacements`` dictionary, leaving other
    types of values (floats, ints, ...) untouched.

    :return: a new list_of_lists with the same shape, with strings replaced by numbers.
    :raise ValueError: if one of the values is not found
    """
    if isinstance(list_of_lists, str):  # Check this first, a str is also Iterable
        try:
            return replacements[list_of_lists]
        except KeyError:
            raise ValueError("Unknown replacement '{}'".format(list_of_lists))
    elif isinstance(list_of_lists, Iterable):
        return [
            replace_symbols_with_values(elem, replacements=replacements)
            for elem in list_of_lists
        ]
    else:
        return list_of_lists  # if it's a numeric value, for instance


def replace_linear_combinations(list_of_3x3_matrices):
    """
    Given a list of 3x3 matrices, where elements can either be float values 
    or lists representing linear combination of values, return a (copied) list
    of 3x3 matrices, where linear combinations are replaced by their values.

    For instance, a value ``[[0.1, 0.3], [0.8, 0.7]]`` means a combination 
    ``0.1 * 0.3 + 0.8 * 0.7`` and will therefore be replaced by ``0.59``.
    """
    result = []

    for matrix in list_of_3x3_matrices:
        new_matrix = []
        for row in matrix:
            new_row = []
            for entry in row:
                if isinstance(entry, Iterable):
                    new_entry = 0
                    for value, factor in entry:
                        new_entry += value * factor
                    new_row.append(new_entry)
                else:
                    new_row.append(entry)
            new_matrix.append(new_row)
        result.append(new_matrix)

    return result
