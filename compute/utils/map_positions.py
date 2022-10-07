import numpy as np
from numpy.random import random


def are_same_except_order(pos1, pos2, cell):
    """
    Pass as pos1 and pos2 the part of the position arrays of two cells (sharing a common cell)
    that have the same species.

    This function returns a tuple ``(distances, pos2_new_order)``, where:
    - ``distances`` is a list (of length N) with the distance between each pair of atoms (after sorting
      the second one according to ``new_order``, in the same Cartesian units as the input)
    - ``pos2_new_order`` a new order for the atoms in pos2, where the first value is the index in pos2 of the atom
       that is closest to the first in pos1, the second value is the index in pos2 of the atom
       that is closest to the second in pos1, ...

    :param pos1: Nx3 array of positions (absolute/Cartesian) for first cell or layer
    :param pos2: Nx3 array of positions (absolute/Cartesian) for second cell or layer
    :param cell: 3x3 array of the lattice vectors (vectors are rows), MUST BE COMMON between the two layers
    """
    assert pos1.shape == pos2.shape, "Error, the two arrays have different shape"

    # Get scaled coordinates
    pos1_scaled = np.dot(pos1, np.linalg.inv(cell))
    pos2_scaled = np.dot(pos2, np.linalg.inv(cell))

    # These will be the output arrays
    distances = []
    pos2_new_order = []

    # This is the array At the beginning I consider all atoms
    remaining = list(range(len(pos2_scaled)))

    # Loop over atoms of pos1, for each find the closest in the "remaining" of pos2
    # (i.e., the ones not already considered in the previous loop iterations)
    for pos in pos1_scaled:
        # Compute difference of scaled positions, modulo 1, making sure that they are in the range [-0.5, 0.5]
        # rather than [0, 1] as it is the default for ``%1``
        pos_diff = (pos - pos2_scaled[np.array(remaining)] + 0.5) % 1.0 - 0.5
        # In theory, only one should be (0, 0, 0). Therefore, I take the minimum row.
        # However, to do this in a Cartesian sense rather than in a scaled sense, I first
        # bring ``pos_diff`` back into Cartesian coordinates
        pos_diff_cart = np.dot(pos_diff, cell)
        # Norm of the distance between two atoms that I believe are the same,
        # after refolding to get the closest
        refolded_norms = (pos_diff_cart**2).sum(axis=1)
        # Note that this index_smallest goes only from 0 to len(remaining)-1
        internal_index_smallest = int(refolded_norms.argmin())
        distances.append(np.sqrt(refolded_norms[internal_index_smallest]))
        # Here I get the correct index in the original array
        index = remaining[internal_index_smallest]
        pos2_new_order.append(index)
        remaining.remove(index)

    # Check that the 'remaining' list is empty.
    # In principle this should never happen, most probably I would get
    # an IndexError in the previous loop
    assert not remaining, "There are some remaining atoms: {}".format(remaining)

    return np.array(distances), np.array(pos2_new_order)


if __name__ == "__main__":
    pos_test = np.array(
        [
            [0.68293989, 0.1153464, 0.9706784],
            [0.62917923, 0.04814019, 0.27771483],
            [0.9845477, 0.81898472, 0.86463285],
            [0.34486299, 0.9378071, 0.49507066],
            [0.21938078, 0.54526453, 0.0271884],
            [0.28520016, 0.37807692, 0.70578199],
            [0.20146576, 0.49949945, 0.72538517],
            [0.75409168, 0.35706838, 0.17894225],
            [0.11575268, 0.18388506, 0.2078017],
            [0.95240967, 0.74938351, 0.57040351],
        ]
    )

    shuffle = np.array([1, 7, 5, 2, 9, 4, 3, 6, 0, 8])

    cell = np.array(
        [
            [3.0, 7.0, 1.0],
            [-5.0, 2.0, 8.0],
            [-2.0, 6.0, -4.0],
        ]
    )

    rel_vec_shift = np.array(
        [
            [1, 1, 1],
            [0, 0, -1],
            [0, 1, 1],
            [0, 1, 1],
            [-1, 0, 0],
            [1, -1, 1],
            [1, 0, 0],
            [-1, 1, 0],
            [-1, -1, -1],
            [-1, 0, -1],
        ]
    )

    NOISE_SIZE = 1.0e-8

    pos_test_compare = (pos_test + rel_vec_shift + NOISE_SIZE * random((10, 3)))[
        shuffle
    ]

    # Add noise?

    pos_test_abs = np.dot(pos_test, cell)
    pos_test_compare_abs = np.dot(pos_test_compare, cell)

    distances, new_order = are_same_except_order(
        pos_test_abs, pos_test_compare_abs, cell
    )
    print("MAX DISTANCE:", distances.max(), "with noise level", NOISE_SIZE)
    assert shuffle[new_order].tolist() == list(range(len(pos_test)))
