import json
import logging
import os
import time
import traceback

import flask
import numpy as np

from .layer_raman_engine import process_structure_core
from .utils.structures import parse_structure
from .utils.matrices import (
    replace_symbols_with_values,
    replace_linear_combinations,
    get_eigvals_eigvects,
)
from .utils.optics import assign_representation, INFRARED, RAMAN, BACKSCATTERING
from .utils.pointgroup import Pointgroup, POINTGROUP_MAPPING
from .utils.response import make_response, FlaskRedirectException


VALID_EXAMPLES = {
    "WTe2": ("WTe2-02f1827d-f339-436f-baf6-66d1cf142fcf_structure.xsf", 1.1),
    "ZnCl2": ("ZnCl2-e5f429a4-3b02-4fb0-8921-0a7ab05078ed_structure.xsf", 1.1),
    # MoS2 bulk from Materials Cloud:
    # https://www.materialscloud.org/explore/2dstructures/details/6e58409f-4ab2-4883-9686-87d4d89c0bf9
    # (Originally from COD, 9007660, P6_3/mmc)
    "MoS2": ("MoS2-6e58409f-4ab2-4883-9686-87d4d89c0bf9_structure.xsf", 1.1),
    # black P bulk from Materials Cloud:
    # https://www.materialscloud.org/explore/2dstructures/details/904c1f0e-da23-42f0-95b4-a4fee98e6d04
    # (Originally from COD, 9012486, Cmce)
    "blackP": ("P-904c1f0e-da23-42f0-95b4-a4fee98e6d04_structure.xsf", 1.1),
    "graphite": ("graphite-544d62e4-8ebe-404c-aa17-b99be62ea70b.xsf", 1.1),
    "BN": ("BN-P6_3mmc-f7e2ff32-27ed-4c89-9c3c-4acbaffbb897.xsf", 1.1),
    "BiTeCl": ("BiTeCl-ICSD79362-category-II.xsf", 1.3),
}


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
            traceback.print_exc()
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
                "user_templates/visualizer.j2", **data_for_template
            )
        except FlaskRedirectException as e:
            flask.flash(str(e))
            return flask.redirect(flask.url_for("input_data"))
        except Exception as exc:
            traceback.print_exc()
            flask.flash(
                "Unable to process the structure, sorry... ({}, {})".format(
                    str(type(exc)), str(exc)
                )
            )
            return flask.redirect(flask.url_for("input_data"))
    else:  # GET Request
        return flask.redirect(flask.url_for("input_data"))


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
                "user_templates/visualizer.j2", **data_for_template
            )
        except FlaskRedirectException as e:
            flask.flash(str(e))
            return flask.redirect(flask.url_for("input_data"))
        except Exception as exc:
            traceback.print_exc()
            flask.flash(
                "Unable to process the structure, sorry... ({}, {})".format(
                    str(type(exc)), str(exc)
                )
            )
            return flask.redirect(flask.url_for("input_data"))
    else:  # GET Request
        return flask.redirect(flask.url_for("input_data"))


@blueprint.route("/api/modes/", methods=["POST"])
def get_modes():  # pylint: disable=too-many-statements
    """Main REST API callback to get the mode frequencies and optical activity."""

    def _get_modes_internal():  # pylint: disable=too-many-locals,too-many-statements,too-many-branches
        """This function implements the actual logic.
        
        It is written as an internal function to allow for profiling.
        """
        start_time = time.time()
        if flask.request.method != "POST":
            return make_response("Only POST method allowed", 405)

        # If not present, assume JSON, even if by default it is usually application/x-www-form-urlencoded
        content_type = flask.request.headers.get("Content-Type", "application/json")
        if content_type.partition(";")[0].strip() != "application/json":
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

        try:
            num_layers_bulk = int(data["numLayersBulk"])
            if num_layers_bulk < 1:
                raise ValueError
        except (KeyError, ValueError):
            return make_response(
                "Invalid request, numLayersBulk value not passed or not in valid range",
                400,
            )

        try:
            pointgroupEven = int(data["pointgroupEven"])
            if pointgroupEven < 1 or pointgroupEven > 32:
                raise KeyError
        except (KeyError, ValueError):
            return make_response(
                "Invalid request, pointgroupEven value not passed or not in valid range",
                400,
            )

        try:
            pointgroupOdd = int(data["pointgroupOdd"])
            if pointgroupOdd < 1 or pointgroupOdd > 32:
                raise KeyError
        except (KeyError, ValueError):
            return make_response(
                "Invalid request, pointgroupOdd value not passed or not in valid range",
                400,
            )

        # at least num_layers_bulk, but also at least 2 layers (otherwise for a single layer we only
        # see the trivial acoustic modes)
        min_num_layers = max(2, num_layers_bulk)

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

        ## LOGIC START ##
        plot_data_x = []
        plot_data_y = []
        is_infrared = []
        is_raman = []
        is_back_scattering = []
        irrep_names = []
        along_x = []
        along_y = []
        along_z = []

        for num_layers in range(min_num_layers, max_layers + 1):
            eigvals, eigvects = get_eigvals_eigvects(num_layers, numeric_matrices)

            pointgroup_number = pointgroupOdd if num_layers % 2 else pointgroupEven
            pointgroup_name = POINTGROUP_MAPPING[pointgroup_number][2]
            pointgroup = Pointgroup(pointgroup_name)

            for eigvec in eigvects.T:  # Skip the first three acoustic modes
                # First axis: layer displacement; second axis: xyz
                # Note: the basis-set order is layer1_x, layer1_y, layer1_z, layer2_x, layer2_y, ...
                displacements = eigvec.reshape((num_layers, 3))

                # TODO: check this assumption!
                # transformation is here always the identity matrix, because
                # we are always rotating the cell so that the stacking axis is along z.
                # In a more general code, it is used to rotate the stacking axis along z.
                irrep_name, activity = assign_representation(
                    displacements, pointgroup, transformation=np.identity(3)
                )
                is_infrared.append(activity[INFRARED])
                is_raman.append(activity[RAMAN])
                # Note: this can be true only when isRamanActive is True, as this means
                # "visible in Raman spectroscopy in a back-scattering geometry".
                # We don't check it, it should be the `assign_representation` function to correctly handle all cases
                is_back_scattering.append(activity[BACKSCATTERING])
                irrep_names.append(irrep_name)
                max_displacement_cart_dir = np.abs(displacements.max(axis=0))
                along_x.append(bool(max_displacement_cart_dir[0] > 1.0e-6))
                along_y.append(bool(max_displacement_cart_dir[1] > 1.0e-6))
                along_z.append(bool(max_displacement_cart_dir[2] > 1.0e-6))

            # 3(N-1) modes, with x equal to the current number of layers (reminder: we are in a loop over num_layers)
            plot_data_x += [num_layers] * 3 * (num_layers - 1)
            plot_data_y += eigvals.tolist()

        return_data = {
            "x": plot_data_x,
            "y": plot_data_y,
            "isBackScattering": is_back_scattering,
            "isRamanActive": is_raman,
            "isInfraredActive": is_infrared,
            "irrepNames": irrep_names,
            "alongX": along_x,
            "alongY": along_y,
            "alongZ": along_z,
        }

        total_time = time.time() - start_time
        return_data["total_time"] = total_time

        ## LOGIC END ##

        print(f"Valid request processed in {total_time} s.")

        return return_data

    ## Uncomment the next four lines (and comment the fifth one) to get also profiling information
    ## This can then be obtained with:
    ## docker cp layer-raman-tool-instance:/home/app/profile.dump .
    ## and then visualized with:
    ## snakeviz profile.dump
    # import cProfile
    # profiler = cProfile.Profile()
    # return_data = profiler.runcall(_get_modes_internal) #, *args, **kwargs)
    # profiler.dump_stats('profile.dump')
    return_data = _get_modes_internal()

    return flask.jsonify(return_data)
