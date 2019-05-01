class FlaskRedirectException(Exception):
    """
    Class used to return immediately with a flash message and a redirect.
    """
    pass

def parse_structure(filecontent, fileformat):
    return (filecontent, fileformat)

def process_structure_core(structure, logger, flask_request, skin_factor):
    """
    Get the data to put in the visualizer jinja2 template
    """


    return {
        'test_data': "Some data from the server-side python code: SLIDER: {} ({}); File format: {}, file content: {} bytes".format(
            skin_factor, str(type(skin_factor)), structure[1], len(structure[0])
        )
    }

