#!/usr/bin/env python
import sys

import ase
import ase.io

from compute.layer_raman_engine import (
    _find_layers,
    find_common_transformation,
    construct_all_matrices,
)


def test_structure(filename, factor=1.1):
    asecell = ase.io.read(filename)
    is_layered, asecell, layer_indices = _find_layers(asecell, factor=factor)
    if not is_layered:
        raise ValueError("The material is not layered")
    rot, transl = find_common_transformation(asecell, layer_indices)
    print("Common transformation found: ")
    print(rot, transl)
    all_matrices = construct_all_matrices(asecell, layer_indices, rot)
    print(all_matrices)


if __name__ == "__main__":
    try:
        filename = sys.argv[1]
    except IndexError:
        print("Pass a filename")
        sys.exit(1)

    try:
        factor = float(sys.argv[2])
    except IndexError:
        factor = 1.1
    except ValueError:
        print("The second parameter, if present, must be a float")
        sys.exit(1)

    test_structure(filename, factor=factor)
