import flask

VERBOSE = True


class FlaskRedirectException(Exception):
    """
    Class used to return immediately with a flash message and a redirect.
    """


def make_response(message, error_code):
    if VERBOSE:
        print(message)
    return flask.make_response(message, error_code)
