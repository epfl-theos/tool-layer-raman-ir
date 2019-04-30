import flask
from flask import Blueprint

import json
import numpy as np

from .layer_raman_engine import process_structure_core, FlaskRedirectException

blueprint = Blueprint('compute', __name__, url_prefix='/compute')

@blueprint.route('/process_structure/', methods=['GET', 'POST'])
def process_structure():
    if flask.request.method == 'POST':
        # check if the post request has the file part
        if 'structurefile' not in flask.request.files:
            return flask.redirect(flask.url_for('input_structure'))
        structurefile = flask.request.files['structurefile']
        fileformat = flask.request.form.get('fileformat', 'unknown')
        filecontent = structurefile.read().decode('utf-8')

        try:
            data_for_template = process_structure_core(
                filecontent=filecontent,
                fileformat=fileformat,
                logger=None,
                flask_request=flask.request)
            return flask.render_template(
                "user_templates/visualizer.html", **data_for_template)
        except FlaskRedirectException as e:
            flask.flash(str(e))
            return flask.redirect(flask.url_for('input_structure'))
        except Exception as exc:
            flask.flash("Unable to process the structure, sorry... ({}, {})".format(str(type(exc)), str(exc)))
            return flask.redirect(flask.url_for('input_data'))
    else:  # GET Request
        return flask.redirect(flask.url_for('input_data'))


verbose = True

def make_response(message, error_code):
    if verbose:
        print(message)
    return flask.make_response(message, error_code)

@blueprint.route('/api/modes/', methods=['POST'])
def get_modes():
    if flask.request.method != 'POST':
        return make_response("Only POST method allowed", 405)
    
    # If not present, assume JSON, even if by default it is usually application/x-www-form-urlencoded
    content_type = flask.request.headers.get('Content-Type', 'application/json') 
    if content_type.partition(';')[0].strip() != 'application/json':    
        print(content_type.partition(';')[0].strip())
        return make_response(
            "Invalid request, Content-Type must be 'application/json' instead of '{}'".format(content_type), 
            400)

    try:
        data = json.loads(flask.request.data.decode('utf-8'))
    except ValueError:
        print(flask.request.data.decode('utf-8'))
        return make_response("Invalid request, not a valid JSON", 400)

    try:
        max_layers = int(data['maxLayers'])
        if max_layers < 2 or max_layers > 40:
            raise KeyError
    except (KeyError, ValueError):
        return make_response("Invalid request, maxLayers value not passed or not in valid range", 400)
        
    try:
        force_constant_params = data['forceConstantParams']
    except KeyError:
        return make_response("Invalid request, missing forceConstantParams", 400)

    # an example of validation
    if 'C111' not in force_constant_params or force_constant_params['C111'] < 0:
        return make_response("missing or invalid C111", 400)

    ## LOGIC START ##

    print(max_layers, force_constant_params)

    num_points = 300
    min_x = -1
    max_x = 1

    x = np.linspace(min_x, max_x, num_points)
    y = np.cos(x * np.pi * max_layers)

    #sleep(0.5)

    return_data = {
            "x": list(x),
            "y": list(y),
            "isBackScattering": [idx % 2 < 1 and (idx % 4 < 2) for idx in range(len(x))], # If it's not Raman active, then this is always False
            "isRamanActive": [idx % 4 < 2 for idx in range(len(x))],
            "isInfraredActive": [idx % 8 < 4 for idx in range(len(x))],
            }

    ## LOGIC END ##

    return flask.jsonify(return_data)