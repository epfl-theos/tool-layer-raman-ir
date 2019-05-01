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
                [[ [['C111', -1.]] , 0., 0.], [0., 'C111', 0.], [0., 0., 'C133']]
            ]
        }
    }

    return {
        'test_data': "Some data from the server-side python code: SLIDER: {}; Number of atoms: {}, chemical numbers: {}".format(
            skin_factor, len(structure[1]), structure[2]
        ),
        'app_data_json': json.dumps(app_data)
    }
