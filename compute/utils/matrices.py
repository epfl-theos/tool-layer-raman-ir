from collections.abc import Iterable
import numpy as np
import scipy
import scipy.linalg


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


def add_block(
    matrix, block, block_i, block_j, factor, banded
):  # pylint: disable=too-many-arguments
    """Add (in place) the block `block` in matrix `matrix`, at block position (i, j) (zero-based).

    The block is first multiplied by `factor` (e.g., a sign).
    If `banded` is True, fill in a banded matrix (with block-band size 1, i.e. the block-diagonal,
    and one upper block-diagonal, and one lower block-diagonal). I.e., `u=1` in the
    documentation of scipy.linalg.eig_banded.
    The banded matrix follows the LAPACK/scipy.linalg.eig_banded structure, when used with 
    `lower=False`. Otherwise, fill in a full matrix.
    
    Matrices must be already initialised with the correct size; the block size is deduced from
    the shape of `block`.
    :param matrix: a numpy array with the matrix to fill
    :param block: a numpy square block
    :para block_i: the row coordinate of the block to add
    :para block_y: the column coordinate of the block to add
    :param factor: a factor to multiply by before adding (factor=1 just adds the block)
    :param banded: if True, fill the position as defined by scipy.linalg.eig_banded, skipping
        elements below the diagonal; otherwise fill a full matrix.
    """
    block_size, block_size_y = block.shape
    assert block_size == block_size_y, "Only square blocks allowed"
    if banded:
        u = block_size
        for i in range(block_size):
            for j in range(block_size):
                actual_i = block_i * block_size + i
                actual_j = block_j * block_size + j
                if actual_i > actual_j:
                    continue
                matrix[u + actual_i - actual_j, actual_j] += factor * block[i, j]
    else:
        matrix[
            block_size * block_i : block_size * (block_i + 1),
            block_size * block_j : block_size * (block_j + 1),
        ] += (factor * block)


def matrix_initialization(variable):
    """
    Initialize the value of a varible inside a 
    force constant matrix, assuming the variable name
    has the form 'Cnab' where n is an integer identifying
    a set of independent variables and a and b identify 
    the directions (1,2 or 3).

    Return a random number between 5 and 10 for the ab=33 component,
    a random number between 1 and 4 for the 11 and 12 components,
    and a random number between -1 and 1 for the off-diagonal components.s
    """
    if variable[-2:] == "33":
        return np.round(5.0 * (1.0 + np.random.rand()), 2)
    if variable[-1] == variable[-2]:
        return np.round(1.0 + 4.0 * np.random.rand(), 2)
    return np.round(-1.0 + 2 * np.random.rand(), 2)


def get_eigvals_eigvects(num_layers, numeric_matrices, use_banded_algorithm=False):
    """Given the number of layers and the `numeric_matrices`,
    construct internally the K matrix and diagonalize it
    to obtain phonon frequencies and modes.

    I have to solve the eq. of motion: omega^2 U = K U

    :param num_layers: the number of layers in the multilayer
    :param numeric_matrices: a list of numeric matrices (i.e., with all parameters replaced with numeric
        values) with the interaction of each layer with the next one.
    :param use_banded_algorithm: if True, use a banded diagonalization algorithm, otherwise diagonalize
        the full matrix. The banded algorithm is always faster and scaled better with the number of
        layers, so there is no reason to set it to False except for debug reasons.

    :return: (eigvals, eigvects), where eigvals is a list of eigenvalues and eigvects a list of list of
        eigenvectors. IMPORTANT! The first three acoustic modes are removed.
        Moreover, the i-th eigenvector is eigenvect.T[i] (note the transpose).
    """
    # TODO: add layer mass here should be a constant for all layers! and fix units

    if use_banded_algorithm:
        # 3 blocks (below, same layer, and above) of size 3 => total width of 9
        # Since we only store the upper part, we only need a width of 4 (diagonal + 3 superdiagonals)
        K_matrix = np.zeros((4, num_layers * 3))
    else:
        K_matrix = np.zeros((num_layers * 3, num_layers * 3))

    for block_idx in range(num_layers):
        # Interaction with upper layer
        if block_idx < num_layers - 1:  # Not in the last layer
            current_block = np.array(
                numeric_matrices[block_idx % len(numeric_matrices)]
            )
            add_block(
                matrix=K_matrix,
                block=current_block,
                block_i=block_idx,
                block_j=block_idx,
                factor=+1,
                banded=use_banded_algorithm,
            )
            add_block(
                matrix=K_matrix,
                block=current_block,
                block_i=block_idx + 1,
                block_j=block_idx,
                factor=-1,
                banded=use_banded_algorithm,
            )
        # Interaction with lower layer
        if block_idx > 0:  # Not in the first layer
            previous_block = np.array(
                numeric_matrices[(block_idx - 1) % len(numeric_matrices)]
            )
            add_block(
                matrix=K_matrix,
                block=previous_block,
                block_i=block_idx,
                block_j=block_idx,
                factor=+1,
                banded=use_banded_algorithm,
            )
            add_block(
                matrix=K_matrix,
                block=previous_block,
                block_i=block_idx - 1,
                block_j=block_idx,
                factor=-1,
                banded=use_banded_algorithm,
            )

    # Get frequencies (eigvals) and eigenvectors (for mode analysis)
    if use_banded_algorithm:
        eigvals, eigvects = scipy.linalg.eig_banded(K_matrix, lower=False)
    else:
        eigvals, eigvects = np.linalg.eigh(K_matrix)

    # The first three should be acoustic i.e. almost zero; the rest should be positive
    assert np.sum(np.abs(eigvals[:3])) < 1.0e-8

    # Remove the first three acoustic modes
    return eigvals[3:], eigvects[:, 3:]
