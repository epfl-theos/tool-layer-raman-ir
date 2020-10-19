from collections.abc import Iterable
import numpy as np

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

def matrix_initialization(variable):
    """
    Initialize the value of a varible inside a 
    force constant matrix, assuming the variable name
    has the form 'Cnab' where n is an integer identifying
    a set of independent variables and a and b identify 
    the directions (1,2 or 3)
    """
    if variable[-2:] == "33":
        return np.round(5.0 * (1.0 + np.random.rand()), 2)
    if variable[-1] == variable[-2]:
        return  np.round(1.0 + 4.0 * np.random.rand(), 2)
    return   np.round(-1.0 + 2 * np.random.rand(), 2)