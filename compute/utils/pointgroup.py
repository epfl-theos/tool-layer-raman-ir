"""
This modules contains all the information concerning pointgroups
that we need to obtain the Raman/Infrared activity of layer modes
"""
import numpy as n
from ase.quaternions import Quaternion
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import copy


SYMPREC = 5.0e-2

# international number, Hermann-Mauguin notation, Schoenflies notation
POINTGROUP_MAPPING = [
    (0, "", ""),  # We put this so we can directly use indexing
    (1, "1", "C1"),
    (2, "-1", "Ci"),
    (3, "2", "C2"),
    (4, "m", "Cs"),
    (5, "2/m", "C2h"),
    (6, "222", "D2"),
    (7, "mm2", "C2v"),
    (8, "mmm", "D2h"),
    (9, "4", "C4"),
    (10, "-4", "S4"),
    (11, "4/m", "C4h"),
    (12, "422", "D4"),
    (13, "4mm", "C4v"),
    (14, "-42m", "D2d"),
    (15, "4/mmm", "D4h"),
    (16, "3", "C3"),
    (17, "-3", "C3i"),
    (18, "32", "D3"),
    (19, "3m", "C3v"),
    (20, "-3m", "D3d"),
    (21, "6", "C6"),
    (22, "-6", "C3h"),
    (23, "6/m", "C6h"),
    (24, "622", "D6"),
    (25, "6mm", "C6v"),
    (26, "-6m2", "D3h"),
    (27, "6/mmm", "D6h"),
    # Below are cubic, will not be really used in this code but we keep them for completeness
    (28, "23", "T"),
    (29, "m-3", "Th"),
    (30, "432", "O"),
    (31, "-43m", "Td"),
    (32, "m-3m", "Oh"),
]


def pg_number_from_hm_symbol(hm_symbol):
    """Return the international number of a given pointgroup, given its H-M symbol."""
    return [_[1] for _ in POINTGROUP_MAPPING].index(hm_symbol)


def rotation(axis, angle):
    """
    Generate the rotation matrix
    given the rotation axis and angle
    """
    norm = n.linalg.norm(axis)
    tol = 2 * n.finfo(n.float).eps
    q = [n.cos(angle / 2)]
    for a in axis:
        q.append(n.sin(angle / 2.0) * a / norm)
    q = Quaternion(q)
    matrix = q.rotation_matrix()
    matrix[abs(matrix) < tol] = 0.0
    return matrix


def prepare_pointgroup(pg_number):
    """Given a pointgroup number (between 1 and 32 inclusive), return a dictionary to be sent to jinja.

    This includes the number itself, and the name both in Hermann-Mauguin notation and in Schönflies notation,
    nicely formatted in HTML (so, will need the `safe` filter)."""
    return {
        "international_number": pg_number,
        "hm_name": POINTGROUP_MAPPING[pg_number][1]
        .replace("-1", '<span style="text-decoration:overline;">1</span>')
        .replace("-3", '<span style="text-decoration:overline;">3</span>')
        .replace("-4", '<span style="text-decoration:overline;">4</span>')
        .replace("-6", '<span style="text-decoration:overline;">6</span>'),
        "schoenflies_name": "{}<sub>{}</sub>".format(
            POINTGROUP_MAPPING[pg_number][2][0], POINTGROUP_MAPPING[pg_number][2][1:]
        ),
    }


def prepare_spacegroup(spg: SpacegroupAnalyzer):
    """Given a spacegroup pypatgen object, return a dictionary to be sent to jinja.

    This includes the international number, and the name in Hermann-Mauguin notation
    nicely formatted in HTML (so, will need the `safe` filter)."""
    return {
        "international_number": spg.get_space_group_number(),
        "hm_name": spg.get_space_group_symbol()
        .replace("-1", '<span style="text-decoration:overline;">1</span>')
        .replace("-3", '<span style="text-decoration:overline;">3</span>')
        .replace("-4", '<span style="text-decoration:overline;">4</span>')
        .replace("-6", '<span style="text-decoration:overline;">6</span>')
        .replace("_1", "<sub>1</sub>")
        .replace("_2", "<sub>2</sub>")
        .replace("_3", "<sub>3</sub>")
        .replace("_4", "<sub>4</sub>")
        .replace("_5", "<sub>5</sub>"),
    }


def rotoreflection(axis, angle):
    """
    Generate the roto-reflection matrix
    given the rotation axis and angle
    (The transformation is nothing
    but a rotation by angle+pi followed by
    inversion)
    """
    # This is a numpy array, so a minus sign in front is valid
    return -rotation(axis, angle + n.pi)  # pylint: disable=invalid-unary-operand-type


# Shorhand notation for complex factors appearing in the character table
w = n.exp(2j * n.pi / 3.0)
w2 = n.exp(-2j * n.pi / 3.0)

pointgroup_dict = {
    # TRICLINIC
    "C1": {
        "classes": {
            0: [n.identity(3)],
        },  # E
        "character_table": {
            "A": [1],
        },
        "raman": ["A"],
        "infrared": ["A"],
        "HM_symbol": "1",
        "backscattering": [["A"], ["A"], ["A"]],
    },
    "Ci": {
        "classes": {
            0: [n.identity(3)],
            1: [-n.identity(3)],
        },  # E  # I
        "character_table": {
            "Ag": [1, 1],
            "Au": [1, -1],
        },
        "raman": ["Ag"],
        "infrared": ["Au"],
        "HM_symbol": "-1",
        "backscattering": [["Ag"], ["Ag"], ["Ag"]],
    },
    # MONOCLINIC
    "C2": {
        "classes": {
            0: [n.identity(3)],
            1: [rotation([0, 0, 1], n.pi)],
        },  # E  # C2
        "character_table": {
            "A": [1, 1],
            "B": [1, -1],
        },
        "raman": ["A", "B"],
        "infrared": ["A", "B"],
        "HM_symbol": "2",
        "backscattering": [["A", "B"], ["A", "B"], ["A"]],
    },
    "Cs": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotoreflection([0, 0, 1], 0.0)],  # sigma_h
        },
        "character_table": {
            "A'": [1, 1],
            "A''": [1, -1],
        },
        "raman": ["A'", "A''"],
        "infrared": ["A'", "A''"],
        "HM_symbol": "m",
        "backscattering": [["A'", "A''"], ["A'", "A''"], ["A'"]],
    },
    "C2h": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], n.pi)],  # C2
            2: [-n.identity(3)],  # I
            3: [rotoreflection([0, 0, 1], 0.0)],  # sigma_h
        },
        "character_table": {
            "Ag": [1, 1, 1, 1],
            "Bg": [1, -1, 1, -1],
            "Au": [1, 1, -1, -1],
            "Bu": [1, -1, -1, 1],
        },
        "raman": ["Ag", "Bg"],
        "infrared": ["Au", "Bu"],
        "HM_symbol": "2/m",
        "backscattering": [["Ag", "Bg"], ["Ag", "Bg"], ["Ag"]],
    },
    # ORTHORHOMBIC
    "D2": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], n.pi)],  # C2z
            2: [rotation([0, 1, 0], n.pi)],  # C2y
            3: [rotation([1, 0, 0], n.pi)],  # C2x
        },
        "character_table": {
            "A": [1, 1, 1, 1],
            "B1": [1, 1, -1, -1],
            "B2": [1, -1, 1, -1],
            "B3": [1, -1, -1, 1],
        },
        "raman": ["A", "B1", "B2", "B3"],
        "infrared": ["B1", "B2", "B3"],
        "HM_symbol": "222",
        "backscattering": [["A", "B3"], ["A", "B2"], ["A", "B1"]],
    },
    "C2v": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], n.pi)],  # C2
            2: [rotoreflection([0, 1, 0], 0.0)],  # sigma_v
            3: [rotoreflection([1, 0, 0], 0.0)],  # sigma_v
        },
        "character_table": {
            "A1": [1, 1, 1, 1],
            "A2": [1, 1, -1, -1],
            "B1": [1, -1, 1, -1],
            "B2": [1, -1, -1, 1],
        },
        "raman": ["A1", "A2", "B1", "B2"],
        "infrared": ["A1", "B1", "B2"],
        "HM_symbol": "mm2",
        "backscattering": [["A1", "B2"], ["A1", "B1"], ["A1", "A2"]],
    },
    "D2h": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], n.pi)],  # C2
            2: [rotation([0, 1, 0], n.pi)],  # C2y
            3: [rotation([1, 0, 0], n.pi)],  # C2x
            4: [-n.identity(3)],  # I
            5: [rotoreflection([0, 0, 1], 0.0)],  # sigma_h
            6: [rotoreflection([0, 1, 0], 0.0)],  # sigma_v
            7: [rotoreflection([1, 0, 0], 0.0)],  # sigma_v
        },
        "character_table": {
            "Ag": [1, 1, 1, 1, 1, 1, 1, 1],
            "B1g": [1, 1, -1, -1, 1, 1, -1, -1],
            "B2g": [1, -1, 1, -1, 1, -1, 1, -1],
            "B3g": [1, -1, -1, 1, 1, -1, -1, 1],
            "Au": [1, 1, 1, 1, -1, -1, -1, -1],
            "B1u": [1, 1, -1, -1, -1, -1, 1, 1],
            "B2u": [1, -1, 1, -1, -1, 1, -1, 1],
            "B3u": [1, -1, -1, 1, -1, 1, 1, -1],
        },
        "raman": ["Ag", "B1g", "B2g", "B3g"],
        "infrared": ["B1u", "B2u", "B3u"],
        "HM_symbol": "mmm",
        "backscattering": [["Ag", "B3g"], ["Ag", "B2g"], ["Ag", "B1g"]],
    },
    # TETRAGONAL
    "C4": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], n.pi)],  # C2
            2: [rotation([0, 0, 1], n.pi / 2.0)],  # C4+
            3: [rotation([0, 0, 1], -n.pi / 2.0)],  # C4-
        },
        "character_table": {
            "A": [1, 1, 1, 1],
            "B": [1, 1, -1, -1],
            "1E": [1, -1, -1j, 1j],
            "2E": [1, -1, 1j, -1j],
        },
        "raman": ["A", "B", "1E", "2E"],
        "infrared": ["A", "1E", "2E"],
        "HM_symbol": "4",
        "backscattering": [["A", "1E", "2E", "B"], ["A", "1E", "2E", "B"], ["A", "B"]],
    },
    "S4": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], n.pi)],  # C2
            2: [rotoreflection([0, 0, 1], n.pi / 2.0)],  # S4+
            3: [rotoreflection([0, 0, 1], -n.pi / 2.0)],  # S4-
        },
        "character_table": {
            "A": [1, 1, 1, 1],
            "B": [1, 1, -1, -1],
            "1E": [1, -1, -1j, 1j],
            "2E": [1, -1, 1j, -1j],
        },
        "raman": ["A", "B", "1E", "2E"],
        "infrared": ["B", "1E", "2E"],
        "HM_symbol": "-4",
        "backscattering": [["A", "1E", "2E", "B"], ["A", "1E", "2E", "B"], ["A", "B"]],
    },
    "C4h": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], n.pi)],  # C2
            2: [rotation([0, 0, 1], n.pi / 2.0)],  # C4+
            3: [rotation([0, 0, 1], -n.pi / 2.0)],  # C4-
            4: [-n.identity(3)],  # I
            5: [rotoreflection([0, 0, 1], 0.0)],  # sigma_h
            6: [rotoreflection([0, 0, 1], -n.pi / 2.0)],  # S4-
            7: [rotoreflection([0, 0, 1], n.pi / 2.0)],  # S4+
        },
        "character_table": {
            "Ag": [1, 1, 1, 1, 1, 1, 1, 1],
            "Bg": [1, 1, -1, -1, 1, 1, -1, -1],
            "1Eg": [1, -1, -1j, 1j, 1, -1, -1j, 1j],
            "2Eg": [1, -1, 1j, -1j, 1, -1, 1j, -1j],
            "Au": [1, 1, 1, 1, -1, -1, -1, -1],
            "Bu": [1, 1, -1, -1, -1, -1, 1, 1],
            "1Eu": [1, -1, -1j, 1j, -1, 1, 1j, -1j],
            "2Eu": [1, -1, 1j, -1j, -1, 1, -1j, 1j],
        },
        "raman": ["Ag", "Bg", "1Eg", "2Eg"],
        "infrared": ["Au", "1Eu", "2Eu"],
        "HM_symbol": "4/m",
        "backscattering": [
            ["Ag", "1Eg", "2Eg", "Bg"],
            ["Ag", "1Eg", "2Eg", "Bg"],
            ["Ag", "Bg"],
        ],
    },
    "D4": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], n.pi)],  # C2
            2: [
                rotation([0, 0, 1], n.pi / 2.0),
                rotation([0, 0, 1], -n.pi / 2.0),
            ],  # C4
            3: [rotation([1, 0, 0], n.pi), rotation([0, 1, 0], n.pi)],  # C2h
            4: [rotation([1, 1, 0], n.pi), rotation([1, -1, 0], n.pi)],  # C2h'
        },
        "character_table": {
            "A1": [1, 1, 1, 1, 1],
            "A2": [1, 1, 1, -1, -1],
            "B1": [1, 1, -1, 1, -1],
            "B2": [1, 1, -1, -1, 1],
            "E": [2, -2, 0, 0, 0],
        },
        "raman": ["A1", "B1", "B2", "E"],
        "infrared": ["A2", "E"],
        "HM_symbol": "422",
        "backscattering": [["A1", "E", "B1"], ["A1", "E", "B1"], ["A1", "B1", "B2"]],
    },
    "C4v": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], n.pi)],  # C2
            2: [
                rotation([0, 0, 1], n.pi / 2.0),
                rotation([0, 0, 1], -n.pi / 2.0),
            ],  # C4
            3: [
                rotoreflection([1, 0, 0], 0.0),
                rotoreflection([0, 1, 0], 0.0),
            ],  # sigma_v
            4: [
                rotoreflection([1, 1, 0], 0.0),
                rotoreflection([1, -1, 0], 0.0),
            ],  # sigma_v'
        },
        "character_table": {
            "A1": [1, 1, 1, 1, 1],
            "A2": [1, 1, 1, -1, -1],
            "B1": [1, 1, -1, 1, -1],
            "B2": [1, 1, -1, -1, 1],
            "E": [2, -2, 0, 0, 0],
        },
        "raman": ["A1", "B1", "B2", "E"],
        "infrared": ["A1", "E"],
        "HM_symbol": "4mm",
        "backscattering": [["A1", "E", "B1"], ["A1", "E", "B1"], ["A1", "B1", "B2"]],
    },
    "D2d": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], n.pi)],  # C2
            2: [
                rotoreflection([0, 0, 1], n.pi / 2.0),
                rotoreflection([0, 0, 1], -n.pi / 2.0),
            ],  # S4
            3: [rotation([1, 0, 0], n.pi), rotation([0, 1, 0], n.pi)],  # C2h
            4: [
                rotoreflection([1, 1, 0], 0.0),
                rotoreflection([1, -1, 0], 0.0),
            ],  # sigma_v'
        },
        "character_table": {
            "A1": [1, 1, 1, 1, 1],
            "A2": [1, 1, 1, -1, -1],
            "B1": [1, 1, -1, 1, -1],
            "B2": [1, 1, -1, -1, 1],
            "E": [2, -2, 0, 0, 0],
        },
        "raman": ["A1", "B1", "B2", "E"],
        "infrared": ["B2", "E"],
        "HM_symbol": "-42m",
        "backscattering": [["A1", "E", "B1"], ["A1", "E", "B1"], ["A1", "B1", "B2"]],
    },
    "D4h": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], n.pi)],  # C2
            2: [
                rotation([0, 0, 1], n.pi / 2.0),
                rotation([0, 0, 1], -n.pi / 2.0),
            ],  # C4
            3: [rotation([1, 0, 0], n.pi), rotation([0, 1, 0], n.pi)],  # C2h
            4: [rotation([1, 1, 0], n.pi), rotation([1, -1, 0], n.pi)],  # C2h'
            5: [-n.identity(3)],  # I
            6: [rotoreflection([0, 0, 1], 0.0)],  # sigma_h
            7: [
                rotoreflection([0, 0, 1], n.pi / 2.0),
                rotoreflection([0, 0, 1], -n.pi / 2.0),
            ],  # S4
            8: [
                rotoreflection([1, 0, 0], 0.0),
                rotoreflection([0, 1, 0], 0.0),
            ],  # sigma_v
            9: [
                rotoreflection([1, 1, 0], 0.0),
                rotoreflection([1, -1, 0], 0.0),
            ],  # sigma_v'
        },
        "character_table": {
            "A1g": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            "A2g": [1, 1, 1, -1, -1, 1, 1, 1, -1, -1],
            "B1g": [1, 1, -1, 1, -1, 1, 1, -1, 1, -1],
            "B2g": [1, 1, -1, -1, 1, 1, 1, -1, -1, 1],
            "Eg": [2, -2, 0, 0, 0, 2, -2, 0, 0, 0],
            "A1u": [1, 1, 1, 1, 1, -1, -1, -1, -1, -1],
            "A2u": [1, 1, 1, -1, -1, -1, -1, -1, 1, 1],
            "B1u": [1, 1, -1, 1, -1, -1, -1, 1, -1, 1],
            "B2u": [1, 1, -1, -1, 1, -1, -1, 1, 1, -1],
            "Eu": [2, -2, 0, 0, 0, -2, 2, 0, 0, 0],
        },
        "raman": ["A1g", "B1g", "B2g", "Eg"],
        "infrared": ["A2u", "Eu"],
        "HM_symbol": "4/mmm",
        "backscattering": [
            ["A1g", "Eg", "B1g"],
            ["A1g", "Eg", "B1g"],
            ["A1g", "B1g", "B2g"],
        ],
    },
    # TRIGONAL
    "C3": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], 2 * n.pi / 3.0)],  # C3+
            2: [rotation([0, 0, 1], -2 * n.pi / 3.0)],  # C3-
        },
        "character_table": {
            "A": [1, 1, 1],
            "1E": [1, w2, w],
            "2E": [1, w, w2],
        },
        "raman": ["A", "1E", "2E"],
        "infrared": ["A", "1E", "2E"],
        "HM_symbol": "3",
        "backscattering": [["A", "1E", "2E"], ["A", "1E", "2E"], ["A", "1E", "2E"]],
    },
    "C3i": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], 2 * n.pi / 3.0)],  # C3+
            2: [rotation([0, 0, 1], -2 * n.pi / 3.0)],  # C3-
            3: [-n.identity(3)],  # I
            4: [rotoreflection([0, 0, 1], -n.pi / 3.0)],  # S6-
            5: [rotoreflection([0, 0, 1], n.pi / 3.0)],  # S6+
        },
        "character_table": {
            "Ag": [1, 1, 1, 1, 1, 1],
            "1Eg": [1, w2, w, 1, w2, w],
            "2Eg": [1, w, w2, 1, w, w2],
            "Au": [1, 1, 1, -1, -1, -1],
            "1Eu": [1, w2, w, -1, -w2, -w],
            "2Eu": [1, w, w2, -1, -w, -w2],
        },
        "raman": ["Ag", "1Eg", "2Eg"],
        "infrared": ["Au", "1Eu", "2Eu"],
        "HM_symbol": "-3",
        "backscattering": [
            ["Ag", "1Eg", "2Eg"],
            ["Ag", "1Eg", "2Eg"],
            ["Ag", "1Eg", "2Eg"],
        ],
    },
    "D3": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [
                rotation([0, 0, 1], 2 * n.pi / 3.0),
                rotation([0, 0, 1], -2 * n.pi / 3.0),
            ],  # C3
            2: [
                rotation([1, 0, 0], n.pi),
                rotation([-1, n.sqrt(3.0), 0], n.pi),
                rotation([-1, -n.sqrt(3.0), 0], n.pi),
            ],  # C2
        },
        "character_table": {
            "A1": [1, 1, 1],
            "A2": [1, 1, -1],
            "E": [2, -1, 0],
        },
        "raman": ["A1", "E"],
        "infrared": ["A2", "E"],
        "HM_symbol": "32",
        "backscattering": [["A1", "E"], ["A1", "E"], ["A1", "E"]],
    },
    "C3v": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [
                rotation([0, 0, 1], 2 * n.pi / 3.0),
                rotation([0, 0, 1], -2 * n.pi / 3.0),
            ],  # C3
            2: [
                rotoreflection([1, 0, 0], 0.0),
                rotoreflection([-1, n.sqrt(3.0), 0], 0.0),
                rotoreflection([-1, -n.sqrt(3.0), 0], 0.0),
            ],  # sigma_v
        },
        "character_table": {
            "A1": [1, 1, 1],
            "A2": [1, 1, -1],
            "E": [2, -1, 0],
        },
        "raman": ["A1", "E"],
        "infrared": ["A1", "E"],
        "HM_symbol": "3m",
        "backscattering": [["A1", "E"], ["A1", "E"], ["A1", "E"]],
    },
    "D3d": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [
                rotation([0, 0, 1], 2 * n.pi / 3.0),
                rotation([0, 0, 1], -2 * n.pi / 3.0),
            ],  # C3
            2: [
                rotation([1, 0, 0], n.pi),
                rotation([-1, n.sqrt(3.0), 0], n.pi),
                rotation([-1, -n.sqrt(3.0), 0], n.pi),
            ],  # C2
            3: [-n.identity(3)],  # I
            4: [
                rotoreflection([0, 0, 1], -n.pi / 3.0),
                rotoreflection([0, 0, 1], n.pi / 3.0),
            ],  # S3
            5: [
                rotoreflection([1, 0, 0], 0.0),
                rotoreflection([-1, n.sqrt(3.0), 0], 0.0),
                rotoreflection([-1, -n.sqrt(3.0), 0], 0.0),
            ],  # sigma_v
        },
        "character_table": {
            "A1g": [1, 1, 1, 1, 1, 1],
            "A2g": [1, 1, -1, 1, 1, -1],
            "Eg": [2, -1, 0, 2, -1, 0],
            "A1u": [1, 1, 1, -1, -1, -1],
            "A2u": [1, 1, -1, -1, -1, 1],
            "Eu": [2, -1, 0, -2, 1, 0],
        },
        "raman": ["A1g", "Eg"],
        "infrared": ["A2u", "Eu"],
        "HM_symbol": "-3m",
        "backscattering": [["A1g", "Eg"], ["A1g", "Eg"], ["A1g", "Eg"]],
    },
    # HEXAGONAL
    "C6": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], n.pi / 3.0)],  # C6+
            2: [rotation([0, 0, 1], 2 * n.pi / 3.0)],  # C3+
            3: [rotation([0, 0, 1], n.pi)],  # C2
            4: [rotation([0, 0, 1], -2 * n.pi / 3.0)],  # C3-
            5: [rotation([0, 0, 1], -n.pi / 3.0)],  # C6-
        },
        "character_table": {
            "A": [1, 1, 1, 1, 1, 1],
            "B": [1, -1, 1, -1, 1, -1],
            "1E2": [1, w, w2, 1, w, w2],
            "2E2": [1, w2, w, 1, w2, w],
            "1E1": [1, -w, w2, -1, w, -w2],
            "2E1": [1, -w2, w, -1, w2, -w],
        },
        "raman": ["A", "1E2", "2E2", "1E1", "2E1"],
        "infrared": ["A", "1E1", "2E1"],
        "HM_symbol": "6",
        "backscattering": [
            ["A", "1E1", "2E1", "1E2", "2E2"],
            ["A", "1E1", "2E1", "1E2", "2E2"],
            ["A", "1E2", "2E2"],
        ],
    },
    "C3h": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], 2 * n.pi / 3.0)],  # C3+
            2: [rotation([0, 0, 1], -2 * n.pi / 3.0)],  # C3-
            3: [rotoreflection([0, 0, 1], 0.0)],  # C2
            4: [rotoreflection([0, 0, 1], -2 * n.pi / 3.0)],  # S3-
            5: [rotoreflection([0, 0, 1], 2 * n.pi / 3.0)],  # S3+
        },
        "character_table": {
            "A'": [1, 1, 1, 1, 1, 1],
            "A''": [1, 1, 1, -1, -1, -1],
            "2E'": [1, w, w2, 1, w, w2],
            "1E'": [1, w2, w, 1, w2, w],
            "2E''": [1, w, w2, -1, -w, -w2],
            "1E''": [1, w2, w, -1, -w2, -w],
        },
        "raman": ["A'", "1E'", "2E'", "1E''", "2E''"],
        "infrared": ["A''", "1E'", "2E'"],
        "HM_symbol": "-6",
        "backscattering": [
            ["A'", "1E''", "2E''", "1E'", "2E'"],
            ["A'", "1E''", "2E''", "1E'", "2E'"],
            ["A'", "1E'", "2E'"],
        ],
    },
    "C6h": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], n.pi / 3.0)],  # C6+
            2: [rotation([0, 0, 1], 2 * n.pi / 3.0)],  # C3+
            3: [rotation([0, 0, 1], n.pi)],  # C2
            4: [rotation([0, 0, 1], -2 * n.pi / 3.0)],  # C3-
            5: [rotation([0, 0, 1], -n.pi / 3.0)],  # C6-
            6: [-n.identity(3)],  # E
            7: [rotoreflection([0, 0, 1], -2 * n.pi / 3.0)],  # S3-
            8: [rotoreflection([0, 0, 1], -n.pi / 3.0)],  # S6-
            9: [rotoreflection([0, 0, 1], 0.0)],  # sigma_h
            10: [rotoreflection([0, 0, 1], n.pi / 3.0)],  # S6+
            11: [rotoreflection([0, 0, 1], 2 * n.pi / 3.0)],  # S3+
        },
        "character_table": {
            "Ag": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            "Bg": [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
            "1E2g": [1, w, w2, 1, w, w2, 1, w, w2, 1, w, w2],
            "2E2g": [1, w2, w, 1, w2, w, 1, w2, w, 1, w2, w],
            "1E1g": [1, -w, w2, -1, w, -w2, 1, -w, w2, -1, w, -w2],
            "2E1g": [1, -w2, w, -1, w2, -w, 1, -w2, w, -1, w2, -w],
            "Au": [1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1],
            "Bu": [1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1],
            "1E2u": [1, w, w2, 1, w, w2, -1, -w, -w2, -1, -w, -w2],
            "2E2u": [1, w2, w, 1, w2, w, -1, -w2, -w, -1, -w2, -w],
            "1E1u": [1, -w, w2, -1, w, -w2, -1, w, -w2, 1, -w, w2],
            "2E1u": [1, -w2, w, -1, w2, -w, -1, w2, -w, 1, -w2, w],
        },
        "raman": ["Ag", "1E2g", "2E2g", "1E1g", "2E1g"],
        "infrared": ["Au", "1E1u", "2E1u"],
        "HM_symbol": "6/m",
        "backscattering": [
            ["Ag", "1E1g", "2E1g", "1E2g", "2E2g"],
            ["Ag", "1E1g", "2E1g", "1E2g", "2E2g"],
            ["Ag", "1E2g", "2E2g"],
        ],
    },
    "D6": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], n.pi)],  # C2
            2: [
                rotation([0, 0, 1], 2 * n.pi / 3.0),
                rotation([0, 0, 1], -2 * n.pi / 3.0),
            ],  # C3
            3: [
                rotation([0, 0, 1], n.pi / 3.0),
                rotation([0, 0, 1], -n.pi / 3.0),
            ],  # C6
            4: [
                rotation([1, 0, 0], n.pi),
                rotation([-1, n.sqrt(3.0), 0], n.pi),
                rotation([-1, -n.sqrt(3.0), 0], n.pi),
            ],  # C2h
            5: [
                rotation([0, 1, 0], n.pi),
                rotation([n.sqrt(3.0), 1, 0], n.pi),
                rotation([n.sqrt(3.0), -1, 0], n.pi),
            ],  # C2h'
        },
        "character_table": {
            "A1": [1, 1, 1, 1, 1, 1],
            "A2": [1, 1, 1, 1, -1, -1],
            "B1": [1, -1, 1, -1, 1, -1],
            "B2": [1, -1, 1, -1, -1, 1],
            "E2": [2, 2, -1, -1, 0, 0],
            "E1": [2, -2, -1, 1, 0, 0],
        },
        "raman": ["A1", "E1", "E2"],
        "infrared": ["A2", "E1"],
        "HM_symbol": "622",
        "backscattering": [["A1", "E1", "E2"], ["A1", "E1", "E2"], ["A1", "E2"]],
    },
    "C6v": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotation([0, 0, 1], n.pi)],  # C2
            2: [
                rotation([0, 0, 1], 2 * n.pi / 3.0),
                rotation([0, 0, 1], -2 * n.pi / 3.0),
            ],  # C3
            3: [
                rotation([0, 0, 1], n.pi / 3.0),
                rotation([0, 0, 1], -n.pi / 3.0),
            ],  # C6
            4: [
                rotoreflection([1, 0, 0], 0.0),
                rotoreflection([-1, n.sqrt(3.0), 0], 0.0),
                rotoreflection([-1, -n.sqrt(3.0), 0], 0.0),
            ],  # sigma_v
            5: [
                rotoreflection([0, 1, 0], 0.0),
                rotoreflection([n.sqrt(3.0), 1, 0], 0.0),
                rotoreflection([n.sqrt(3.0), -1, 0], 0.0),
            ],  # sigma_v'
        },
        "character_table": {
            "A1": [1, 1, 1, 1, 1, 1],
            "A2": [1, 1, 1, 1, -1, -1],
            "B1": [1, -1, 1, -1, 1, -1],
            "B2": [1, -1, 1, -1, -1, 1],
            "E2": [2, 2, -1, -1, 0, 0],
            "E1": [2, -2, -1, 1, 0, 0],
        },
        "raman": ["A1", "E1", "E2"],
        "infrared": ["A1", "E1"],
        "HM_symbol": "6mm",
        "backscattering": [["A1", "E1", "E2"], ["A1", "E1", "E2"], ["A1", "E2"]],
    },
    "D3h": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [rotoreflection([0, 0, 1], 0.0)],  # sigma_h
            2: [
                rotation([0, 0, 1], 2 * n.pi / 3.0),
                rotation([0, 0, 1], -2 * n.pi / 3.0),
            ],  # C3
            3: [
                rotoreflection([0, 0, 1], 2 * n.pi / 3.0),
                rotoreflection([0, 0, 1], -2 * n.pi / 3.0),
            ],  # S3
            4: [
                rotation([1, 0, 0], n.pi),
                rotation([-1, n.sqrt(3.0), 0], n.pi),
                rotation([-1, -n.sqrt(3.0), 0], n.pi),
            ],  # C2h
            5: [
                rotoreflection([0, 1, 0], 0.0),
                rotoreflection([n.sqrt(3.0), 1, 0], 0.0),
                rotoreflection([n.sqrt(3.0), -1, 0], 0.0),
            ],  # sigma_v'
        },
        "character_table": {
            "A1'": [1, 1, 1, 1, 1, 1],
            "A2'": [1, 1, 1, 1, -1, -1],
            "A1''": [1, -1, 1, -1, 1, -1],
            "A2''": [1, -1, 1, -1, -1, 1],
            "E'": [2, 2, -1, -1, 0, 0],
            "E''": [2, -2, -1, 1, 0, 0],
        },
        "raman": ["A1'", "E'", "E''"],
        "infrared": ["A2''", "E'"],
        "HM_symbol": "-6m2",
        "backscattering": [["A1'", "E'", "E''"], ["A1'", "E'", "E''"], ["A1'", "E'"]],
    },
    "D6h": {
        "classes": {
            0: [n.identity(3)],  # E
            1: [
                rotation([0, 0, 1], n.pi / 3.0),
                rotation([0, 0, 1], -n.pi / 3.0),
            ],  # C6
            2: [
                rotation([0, 0, 1], 2 * n.pi / 3.0),
                rotation([0, 0, 1], -2 * n.pi / 3.0),
            ],  # C3
            3: [rotation([0, 0, 1], n.pi)],  # C2
            4: [
                rotation([1, 0, 0], n.pi),
                rotation([-1, n.sqrt(3.0), 0], n.pi),
                rotation([-1, -n.sqrt(3.0), 0], n.pi),
            ],  # C2h
            5: [
                rotation([0, 1, 0], n.pi),
                rotation([n.sqrt(3.0), 1, 0], n.pi),
                rotation([n.sqrt(3.0), -1, 0], n.pi),
            ],  # C2h'
            6: [-n.identity(3)],  # I
            7: [
                rotoreflection([0, 0, 1], 2 * n.pi / 3.0),
                rotoreflection([0, 0, 1], -2 * n.pi / 3.0),
            ],  # S3
            8: [
                rotoreflection([0, 0, 1], n.pi / 3.0),
                rotoreflection([0, 0, 1], -n.pi / 3.0),
            ],  # S6
            9: [rotoreflection([0, 0, 1], 0.0)],  # sigma_h
            10: [
                rotoreflection([1, 0, 0], 0.0),
                rotoreflection([-1, n.sqrt(3.0), 0], 0.0),
                rotoreflection([-1, -n.sqrt(3.0), 0], 0.0),
            ],  # sigma_v
            11: [
                rotoreflection([0, 1, 0], 0.0),
                rotoreflection([n.sqrt(3.0), 1, 0], 0.0),
                rotoreflection([n.sqrt(3.0), -1, 0], 0.0),
            ],  # sigma_v'
        },
        "character_table": {
            "A1g": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            "A2g": [1, 1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1],
            "B1g": [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
            "B2g": [1, -1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1],
            "E2g": [2, -1, -1, 2, 0, 0, 2, -1, -1, 2, 0, 0],
            "E1g": [2, 1, -1, -2, 0, 0, 2, 1, -1, -2, 0, 0],
            "A1u": [1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1],
            "A2u": [1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1],
            "B1u": [1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1],
            "B2u": [1, -1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1],
            "E2u": [2, -1, -1, 2, 0, 0, -2, 1, 1, -2, 0, 0],
            "E1u": [2, 1, -1, -2, 0, 0, -2, -1, 1, 2, 0, 0],
        },
        "raman": ["A1g", "E1g", "E2g"],
        "infrared": ["A2u", "E1u"],
        "HM_symbol": "6/mmm",
        "backscattering": [
            ["A1g", "E1g", "E2g"],
            ["A1g", "E1g", "E2g"],
            ["A1g", "E2g"],
        ],
    },
}


class Pointgroup:
    def __init__(self, pointgroupname):
        info = copy.deepcopy(pointgroup_dict[pointgroupname])
        self.symbol = pointgroupname

        self.classes = info.pop("classes")
        self.character_table = info.pop("character_table")
        self.raman = info.pop("raman")
        self.infrared = info.pop("infrared")
        self.HM_symbol = info.pop("HM_symbol")
        self.backscattering = info.pop("backscattering")
        if info:
            raise ValueError(
                "These keys were unexpected in pointgroup_dict['{}']: {}".format(
                    pointgroupname, ", ".join(info.keys())
                )
            )

    def get_order(self):
        """
        Return the order of the group
        by counting the number of elements in each class
        """
        count = 0
        for _, op in self.classes.items():
            count += len(op)
        return count

    def _check_order_by_irreps(self):
        """
        Check that the order obtained from the dimensions of the irreps
        is the same obtained from the elements of the group.

        Here we assume that the first operation (first value of the characters
        list) is always the identity.
        """
        count = 0
        for (
            irreps,  # pylint: disable=unused-variable
            character,
        ) in self.character_table.items():
            count += character[0] ** 2

        return count == self.get_order()

    def _check_orthogonality_irreps(self):
        """
        Check the orthogonality of the irreps
        """
        num_irreps = len(self.character_table)
        test = n.identity(num_irreps) * (self.get_order() + 0j)
        for c, sym_ops in self.classes.items():
            for i, chars1 in enumerate(self.character_table.values()):
                for j, chars2 in enumerate(self.character_table.values()):
                    test[i, j] -= chars1[c] * n.conj(chars2[c]) * len(sym_ops)
        return n.sum(abs(test)) < 1.0e-5

    def _check_orthogonality_classes(self):
        """
        Check orthogonality of the classes
        """
        num_classes = len(self.classes)
        test = n.diag(
            [self.get_order() * (1.0 + 0j) / len(c) for c in self.classes.values()]
        )
        for (
            irreps,  # pylint: disable=unused-variable
            characters,
        ) in self.character_table.items():
            for i in range(num_classes):
                for j in range(num_classes):
                    test[i, j] -= characters[i] * n.conj(characters[j])
        return n.sum(abs(test)) < 1.0e-5

    def check_all(self, do_print=True):
        if do_print:
            print(self.symbol)
        mapping = {True: "Ok", False: "FAIL!"}
        check_result = self._check_order_by_irreps()
        if do_print:
            print("Check order: {}".format(mapping[check_result]))
        else:
            assert check_result

        check_result = self._check_orthogonality_irreps()
        if do_print:
            print("Check orthogonality (wrt irreps): {}".format(mapping[check_result]))
        else:
            assert check_result

        check_result = self._check_orthogonality_classes()
        if do_print:
            print("Check orthogonality (wrt classes): {}".format(mapping[check_result]))
        else:
            assert check_result
