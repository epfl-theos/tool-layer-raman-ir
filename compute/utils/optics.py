"""This module contains the functions used to check the optical activity of a mode."""
import numpy as n


RAMAN = "raman"
INFRARED = "infrared"
BACKSCATTERING = "backscattering"


def symmetry_operation(matrix, displacements):
    """
    This function specifies how a symmetry operation operates on the layers
    """
    tol = 1.0e-5
    # First, we apply the matrix to the displacement of each layer
    new_displacements = n.dot(matrix, displacements.transpose()).transpose()

    # Then we see how the operation rearranges the order of the layers
    # NOTE: there is only one possible centre for this operation in a finite layer, at the centre of the layer,
    # so 1, 2, 3, 4 becomes 4, 3, 2, 1; it cannot become e.g. 2, 1, 4, 3 as it would also shift the finite ML
    if abs(matrix[2, 2] - 1) < tol:
        return new_displacements
    if abs(matrix[2, 2] + 1) < tol:
        return new_displacements[::-1, :]
    raise ValueError(
        ("Transformation not allowed as it does not leave the layers invariant!")
    )


def assign_representation(
    displacements, pg, transformation, do_print=False
):  # pylint: disable=too-many-locals, too-many-branches
    """ 
    This function assigns a set of layer displacements to a given 
    irreducible representation of the pointgroup pg.

    :param displacements: a num_layer x 3 array of displacements
    :param pg: the point group (an instance of utils.Pointgroup)
    :param transformation: a 3x3 matrix that transforms the axes by swapping them
    """
    tol = 1e-5
    idir = n.argwhere(transformation[2, :] == 1)[0, 0]
    activity = {RAMAN: False, INFRARED: False, BACKSCATTERING: False}
    found = False
    foundtwice = False
    for irrep_idx, (irrep, characters) in enumerate(pg.character_table.items()):
        # Determine if the irreps is complex by looking at the characters
        iscomplex = n.max([abs(n.imag(c)) for c in characters]) > 0
        proj = 0.0
        for c, sym_ops in pg.classes.items():
            for sym_op in sym_ops:
                this_sym_op = n.dot(
                    n.dot(transformation.transpose(), sym_op), transformation
                )
                if do_print and irrep_idx == 0:
                    print("this_sym_op:")
                    print(this_sym_op)
                    print("displacements:")
                    print(displacements)
                    print("-" * 20)
                proj += (
                    # complex representations should appear together
                    # because of time reversal symmetry, so we project on
                    # both partner representations
                    (n.conj(characters[c]) + iscomplex * characters[c])
                    * n.diag(
                        n.dot(
                            displacements,
                            symmetry_operation(this_sym_op, displacements).transpose(),
                        )
                    ).sum()
                )
        # We assume the identity to be always the first class
        proj *= characters[0] * 1.0 / pg.get_order()
        # print("   -Irrep: {}, proj: {}".format(irrep,proj))
        if abs(proj) > tol:
            # complex representations should appear twice
            if found:
                if iscomplex:
                    if foundtwice:
                        print("!!!! Found more than twice complex irrep")
                    foundtwice = True
                else:
                    print("!!!! Found more than once non-complex irrep")
            found = True
            # check that the projection is 1, within a threshold
            if abs(proj - 1) > tol:
                print("!!!! Not one or zero, {}, {}, {}".format(proj, irrep, pg.symbol))
            # To avoid issues in case of accidental degeneracies
            if proj > 0.5:
                this_irrep = irrep
                if do_print:
                    print("      {} irrep,".format(irrep), end="")
                if irrep in pg.raman:
                    activity[RAMAN] = True
                    if do_print:
                        print(" Raman+IR active", end="")
                    if irrep in pg.backscattering[idir]:
                        activity[BACKSCATTERING] = True
                if irrep in pg.infrared:
                    activity[INFRARED] = True
                    if do_print:
                        print(" IR active", end="")
                if do_print:
                    print("")
    return this_irrep, activity
