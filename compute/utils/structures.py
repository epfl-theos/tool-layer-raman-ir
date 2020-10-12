import io

import numpy as np
import ase

from tools_barebone.structure_importers import get_structure_tuple, UnknownFormatError
from .response import FlaskRedirectException


def parse_structure(filecontent, fileformat, extra_data=None):
    """Parse a structure given the file content and the file format.

    Possibly pass also the extra data from the form if needed
    (e.g. for the cell of a XYZ file).
    """
    fileobject = io.StringIO(str(filecontent))
    try:
        structure_tuple = get_structure_tuple(
            fileobject, fileformat, extra_data=extra_data
        )
    except UnknownFormatError:
        raise FlaskRedirectException("Unknown format '{}'".format(fileformat))
    # Bubble up the exception, will be managed by the level above

    return structure_tuple


def ase_from_tuple(structure_tuple):
    """
    Convert a structure tuple with the format 
    (cell, scaled_positions,element_numbers) to a ASE structure
    
    Currently, if the element numbers do not correspond to a known chemical species,
    an exception is raised. 

    :param structure_tuple: the structure in format (cell, scaled_positions,element_numbers)
    """
    cell, rel_pos, numbers = structure_tuple

    if any([i not in ase.atom.atomic_numbers.values() for i in numbers]):
        raise ValueError

    return ase.Atoms(cell=cell, scaled_positions=rel_pos, pbc=True, numbers=numbers)


def tuple_from_ase(asecell: ase.Atoms):
    """
    Convert an ase cell to a structure tuple.
    """
    cell = asecell.cell.tolist()
    # Wrap=False to preserve the absolute positions of atoms
    rel_pos = asecell.get_scaled_positions(wrap=False).tolist()
    numbers = [
        ase.atom.atomic_numbers[symbol] for symbol in asecell.get_chemical_symbols()
    ]
    return (cell, rel_pos, numbers)


def get_xsf_structure(structure_tuple):
    """Return a XSF file from a structure tuple."""
    cell = np.array(structure_tuple[0])

    xsfstructure = []
    xsfstructure.append("CRYSTAL")
    xsfstructure.append("PRIMVEC")
    for vector in cell:
        xsfstructure.append("{} {} {}".format(vector[0], vector[1], vector[2]))
    xsfstructure.append("PRIMCOORD")
    xsfstructure.append("{} 1".format(len(structure_tuple[1])))
    for scaled_pos, atom_num in zip(structure_tuple[1], structure_tuple[2]):
        abs_pos = cell.T @ scaled_pos
        xsfstructure.append(
            "{} {} {} {}".format(atom_num, abs_pos[0], abs_pos[1], abs_pos[2])
        )
    xsfstructure = "\n".join(xsfstructure)

    return xsfstructure


def get_covalent_radii_array(asecell):
    """ 
    Return a list of the covalent radii for 
    the atoms in the ASE structure using the Cordero values

    :params asecell: the ASE structure
    """
    map_atomic_number_covalent_cordero = {
        1: 0.31,
        2: 0.28,
        3: 1.28,
        4: 0.96,
        5: 0.84,
        6: 0.76,
        7: 0.71,
        8: 0.66,
        9: 0.57,
        10: 0.58,
        11: 1.66,
        12: 1.41,
        13: 1.21,
        14: 1.11,
        15: 1.07,
        16: 1.05,
        17: 1.02,
        18: 1.06,
        19: 2.03,
        20: 1.76,
        21: 1.7,
        22: 1.6,
        23: 1.53,
        24: 1.39,
        25: 1.39,
        26: 1.32,
        27: 1.26,
        28: 1.24,
        29: 1.32,
        30: 1.22,
        31: 1.22,
        32: 1.2,
        33: 1.19,
        34: 1.2,
        35: 1.2,
        36: 1.16,
        37: 2.2,
        38: 1.95,
        39: 1.9,
        40: 1.75,
        41: 1.64,
        42: 1.54,
        43: 1.47,
        44: 1.46,
        45: 1.42,
        46: 1.39,
        47: 1.45,
        48: 1.44,
        49: 1.42,
        50: 1.39,
        51: 1.39,
        52: 1.38,
        53: 1.39,
        54: 1.4,
        55: 2.44,
        56: 2.15,
        57: 2.07,
        58: 2.04,
        59: 2.03,
        60: 2.01,
        61: 1.99,
        62: 1.98,
        63: 1.98,
        64: 1.96,
        65: 1.94,
        66: 1.92,
        67: 1.92,
        68: 1.89,
        69: 1.9,
        70: 1.87,
        71: 1.87,
        72: 1.75,
        73: 1.7,
        74: 1.62,
        75: 1.51,
        76: 1.44,
        77: 1.41,
        78: 1.36,
        79: 1.36,
        80: 1.32,
        81: 1.45,
        82: 1.46,
        83: 1.48,
        84: 1.4,
        85: 1.5,
        86: 1.5,
        87: 2.6,
        88: 2.21,
        89: 2.15,
        90: 2.06,
        91: 2,
        92: 1.96,
        93: 1.9,
        94: 1.87,
        95: 1.8,
        96: 1.69,
    }
    return np.array(
        [map_atomic_number_covalent_cordero.get(atom.number) for atom in asecell]
    )
