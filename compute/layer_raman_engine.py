class FlaskRedirectException(Exception):
    """
    Class used to return immediately with a flash message and a redirect.
    """
    pass

def process_structure_core(filecontent, fileformat, logger, flask_request, example_slider):
    """
    Get the data to put in the visualizer jinja2 template
    """
    return {
        'test_data': "Some data from the server-side python code: SLIDER: {}; File format: {}, file content: {} bytes".format(
            example_slider, fileformat, len(filecontent)
        )
    }
