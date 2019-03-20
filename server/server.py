#!/usr/bin/env python
from __future__ import print_function
from flask import Flask
import flask
import json
import numpy as np
from flask import request
from flask_cors import CORS

app = Flask(__name__)
CORS(app, support_credentials=True)

verbose = True

def make_response(message, error_code):
    if verbose:
        print(message)
    return flask.make_response(message, error_code)

@app.route('/api/modes/', methods=['POST'])
def get_modes():
    if request.method != 'POST':
        return make_response("Only POST method allowed", 405)
    
    # If not present, assume JSON, even if by default it is usually application/x-www-form-urlencoded
    content_type = request.headers.get('Content-Type', 'application/json') 
    if content_type.partition(';')[0].strip() != 'application/json':    
        print(content_type.partition(';')[0].strip())
        return make_response(
            "Invalid request, Content-Type must be 'application/json' instead of '{}'".format(content_type), 
            400)

    try:
        data = json.loads(request.data.decode('utf-8'))
    except ValueError:
        print(request.data.decode('utf-8'))
        return make_response("Invalid request, not a valid JSON", 400)

    try:
        max_layers = data['maxLayers']
    except KeyError:
        return make_response("Invalid request, maxLayers value not passed", 400)
        
    try:
        max_layers = int(max_layers)
        if max_layers < 2 or max_layers > 40:
            raise KeyError
    except KeyError:
        return make_response("Invalid request, k is not an int value or it is not in the expected range", 400)

    ## LOGIC START ##

    num_points = 300
    min_x = -1
    max_x = 1

    x = np.linspace(min_x, max_x, num_points)
    y = np.cos(x * np.pi * max_layers)

    #sleep(0.5)

    return_data = {
            "x": list(x),
            "y": list(y),
            "isBackScattering": [True for _ in x],
            "isRamanActive": [True for _ in x],
            "isInfraredActive": [True for _ in x],
            }

    ## LOGIC END ##

    return flask.jsonify(return_data)

if __name__ == '__main__':
    # Must be inside this 'if' if you want to use apache mod_wsgi
    # To install mod_wsgi, follow instructions here
    # http://flask.pocoo.org/docs/0.10/deploying/mod_wsgi/
    app.run(port=5300)
