import numpy as np
from numpy import ndarray
from typing import List, Callable, Dict


def get_number_of_quadrature_points_in_segment(integration_order: int) -> int:
    """

    Args:
        integration_order:

    Returns:

    """
    if integration_order in [0, 1]:
        number_of_quadrature_points = 1
    elif integration_order == 2:
        number_of_quadrature_points = 2
    elif integration_order == 3:
        number_of_quadrature_points = 3
    elif integration_order == 4:
        number_of_quadrature_points = 4
    elif integration_order == 5:
        number_of_quadrature_points = 5
    elif integration_order == 6:
        number_of_quadrature_points = 6
    elif integration_order == 7:
        number_of_quadrature_points = 7
    elif integration_order == 8:
        number_of_quadrature_points = 8
    else:
        raise ValueError("quadrature order not supported")
    return number_of_quadrature_points


def get_segment_quadrature_data(integration_order: int,) -> (ndarray, ndarray):
    """

    Args:
        integration_order:

    Returns:

    """
    if integration_order in [0, 1]:
        reference_points = [[+0.5000000000, +0.5000000000]]
        reference_points = [[0.5000000000000000, 0.5000000000000000]]
        # reference_weights = [+0.5000000000]
        reference_weights = [+1.0000000000]
    elif integration_order == 2:
        reference_points = [
            [+0.7886751346, +0.2113248654],
            [+0.2113248654, +0.7886751346],
        ]
        reference_points = [
            [0.78867513459481290, 0.21132486540518713],
            [0.21132486540518713, 0.78867513459481290],
        ]
        reference_weights = [+0.5000000000, +0.5000000000]
    elif integration_order == 3:
        reference_points = [
            [+0.8872983346, +0.1127016654],
            [+0.5000000000, +0.5000000000],
            [+0.1127016654, +0.8872983346],
        ]
        reference_weights = [+0.2777777778, +0.4444444444, +0.2777777778]
    elif integration_order == 4:
        reference_points = [
            [+0.9305681558, +0.0694318442],
            [+0.6699905218, +0.3300094782],
            [+0.3300094782, +0.6699905218],
            [+0.0694318442, +0.9305681558],
        ]
        reference_weights = [+0.1739274226, +0.3260725774, +0.3260725774, +0.1739274226]
    elif integration_order == 5:
        reference_points = [
            [+0.9530899230, +0.0469100770],
            [+0.7692346551, +0.2307653449],
            [+0.5000000000, +0.5000000000],
            [+0.2307653449, +0.7692346551],
            [+0.0469100770, +0.9530899230],
        ]
        reference_weights = [
            +0.1184634425,
            +0.2393143352,
            +0.2844444444,
            +0.2393143352,
            +0.1184634425,
        ]
    elif integration_order == 6:
        reference_points = [
            [+0.9662347571, +0.0337652429],
            [+0.8306046932, +0.1693953068],
            [+0.6193095930, +0.3806904070],
            [+0.3806904070, +0.6193095930],
            [+0.1693953068, +0.8306046932],
            [+0.0337652429, +0.9662347571],
        ]
        reference_weights = [
            +0.0856622462,
            +0.1803807865,
            +0.2339569673,
            +0.2339569673,
            +0.1803807865,
            +0.0856622462,
        ]
    elif integration_order == 7:
        reference_points = [
            [+0.9745539562, +0.0254460438],
            [+0.8707655928, +0.1292344072],
            [+0.7029225757, +0.2970774243],
            [+0.5000000000, +0.5000000000],
            [+0.2970774243, +0.7029225757],
            [+0.1292344072, +0.8707655928],
            [+0.0254460438, +0.9745539562],
        ]
        reference_weights = [
            +0.0647424831,
            +0.1398526957,
            +0.1909150253,
            +0.2089795918,
            +0.1909150253,
            +0.1398526957,
            +0.0647424831,
        ]
    elif integration_order == 8:
        reference_points = [
            [+0.9801449282, +0.0198550718],
            [+0.8983332387, +0.1016667613],
            [+0.7627662050, +0.2372337950],
            [+0.5917173212, +0.4082826788],
            [+0.4082826788, +0.5917173212],
            [+0.2372337950, +0.7627662050],
            [+0.1016667613, +0.8983332387],
            [+0.0198550718, +0.9801449282],
        ]
        reference_weights = [
            +0.0506142681,
            +0.1111905172,
            +0.1568533229,
            +0.1813418917,
            +0.1813418917,
            +0.1568533229,
            +0.1111905172,
            +0.0506142681,
        ]
    else:
        raise ValueError("quadrature order not supported")
    return np.array(reference_points), np.array(reference_weights)
