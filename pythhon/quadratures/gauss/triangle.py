from numpy import ndarray

from pythhon.parameters import *


def get_number_of_quadrature_points_in_triangle(integration_order: int) -> int:
    """

    Args:
        integration_order:

    Returns:

    """
    if integration_order in [0, 1]:
        number_of_quadrature_points = 1
    elif integration_order == 2:
        number_of_quadrature_points = 3
    elif integration_order == 3:
        number_of_quadrature_points = 4
    elif integration_order == 4:
        number_of_quadrature_points = 6
    elif integration_order == 5:
        number_of_quadrature_points = 7
    elif integration_order == 6:
        number_of_quadrature_points = 12
    elif integration_order == 7:
        number_of_quadrature_points = 13
    elif integration_order == 8:
        number_of_quadrature_points = 16
    else:
        raise ValueError("quadrature order not supported")
    return number_of_quadrature_points


def get_triangle_quadrature_data(integration_order: int,) -> (ndarray, ndarray):
    """

    Args:
        integration_order:

    Returns:

    """
    if integration_order in [0, 1]:
        reference_points = [[0.3333333333333333, 0.3333333333333333, 0.3333333333333333]]
        reference_weights = [+1.0000000000]
    elif integration_order == 2:
        # reference_points = [
        #     [+0.6666666667, +0.1666666667, +0.1666666667],
        #     [+0.1666666667, +0.6666666667, +0.1666666667],
        #     [+0.1666666667, +0.1666666667, +0.6666666667],
        # ]
        reference_points = [
            [0.66666666666666666, 0.16666666666666666, 0.16666666666666666],
            [0.16666666666666666, 0.66666666666666666, 0.16666666666666666],
            [0.16666666666666666, 0.16666666666666666, 0.66666666666666666],
        ]
        reference_weights = [+0.3333333333, +0.3333333333, +0.3333333333]
    elif integration_order == 3:
        reference_points = [
            [+0.3333333333, +0.3333333333, +0.3333333333],
            [+0.2000000000, +0.6000000000, +0.2000000000],
            [+0.2000000000, +0.2000000000, +0.6000000000],
            [+0.6000000000, +0.2000000000, +0.2000000000],
        ]
        reference_weights = [-0.5625000000, +0.5208333333, +0.5208333333, +0.5208333333]
    elif integration_order == 4:
        reference_points = [
            [+0.0915762135, +0.8168475730, +0.0915762135],
            [+0.0915762135, +0.0915762135, +0.8168475730],
            [+0.8168475730, +0.0915762135, +0.0915762135],
            [+0.4459484909, +0.1081030182, +0.4459484909],
            [+0.4459484909, +0.4459484909, +0.1081030182],
            [+0.1081030182, +0.4459484909, +0.4459484909],
        ]
        reference_weights = [
            +0.1099517437,
            +0.1099517437,
            +0.1099517437,
            +0.2233815897,
            +0.2233815897,
            +0.2233815897,
        ]
    elif integration_order == 5:
        reference_points = [
            [+0.3333333333, +0.3333333333, +0.3333333333],
            [+0.1012865073, +0.7974269854, +0.1012865073],
            [+0.1012865073, +0.1012865073, +0.7974269854],
            [+0.7974269854, +0.1012865073, +0.1012865073],
            [+0.4701420641, +0.4701420641, +0.0597158718],
            [+0.4701420641, +0.0597158718, +0.4701420641],
            [+0.0597158718, +0.4701420641, +0.4701420641],
        ]
        reference_weights = [
            +0.2250000000,
            +0.1259391805,
            +0.1259391805,
            +0.1259391805,
            +0.1323941528,
            +0.1323941528,
            +0.1323941528,
        ]
    elif integration_order == 6:
        reference_points = [
            [+0.0630890145, +0.8738219710, +0.0630890145],
            [+0.0630890145, +0.0630890145, +0.8738219710],
            [+0.8738219710, +0.0630890145, +0.0630890145],
            [+0.2492867452, +0.5014265097, +0.2492867452],
            [+0.2492867452, +0.2492867452, +0.5014265097],
            [+0.5014265097, +0.2492867452, +0.2492867452],
            [+0.0531450498, +0.6365024991, +0.3103524510],
            [+0.0531450498, +0.3103524510, +0.6365024991],
            [+0.3103524510, +0.6365024991, +0.0531450498],
            [+0.6365024991, +0.3103524510, +0.0531450498],
            [+0.6365024991, +0.0531450498, +0.3103524510],
            [+0.3103524510, +0.0531450498, +0.6365024991],
        ]
        reference_weights = [
            +0.0508449064,
            +0.0508449064,
            +0.0508449064,
            +0.1167862757,
            +0.1167862757,
            +0.1167862757,
            +0.0828510756,
            +0.0828510756,
            +0.0828510756,
            +0.0828510756,
            +0.0828510756,
            +0.0828510756,
        ]
    elif integration_order == 7:
        reference_points = [
            [+0.3333333333, +0.3333333333, +0.3333333333],
            [+0.2603459661, +0.4793080678, +0.2603459661],
            [+0.2603459661, +0.2603459661, +0.4793080678],
            [+0.4793080678, +0.2603459661, +0.2603459661],
            [+0.0651301029, +0.8697397942, +0.0651301029],
            [+0.0651301029, +0.0651301029, +0.8697397942],
            [+0.8697397942, +0.0651301029, +0.0651301029],
            [+0.6384441886, +0.0486903154, +0.3128654960],
            [+0.6384441886, +0.3128654960, +0.0486903154],
            [+0.3128654960, +0.6384441886, +0.0486903154],
            [+0.3128654960, +0.0486903154, +0.6384441886],
            [+0.0486903154, +0.3128654960, +0.6384441886],
            [+0.0486903154, +0.6384441886, +0.3128654960],
        ]
        reference_weights = [
            -0.1495700445,
            +0.1756152574,
            +0.1756152574,
            +0.1756152574,
            +0.0533472356,
            +0.0533472356,
            +0.0533472356,
            +0.0771137609,
            +0.0771137609,
            +0.0771137609,
            +0.0771137609,
            +0.0771137609,
            +0.0771137609,
        ]
    elif integration_order == 8:
        reference_points = [
            [+0.3333333333, +0.3333333333, +0.3333333333],
            [+0.4592925883, +0.0814148234, +0.4592925883],
            [+0.4592925883, +0.4592925883, +0.0814148234],
            [+0.0814148234, +0.4592925883, +0.4592925883],
            [+0.1705693078, +0.6588613845, +0.1705693078],
            [+0.1705693078, +0.1705693078, +0.6588613845],
            [+0.6588613845, +0.1705693078, +0.1705693078],
            [+0.0505472283, +0.8989055434, +0.0505472283],
            [+0.0505472283, +0.0505472283, +0.8989055434],
            [+0.8989055434, +0.0505472283, +0.0505472283],
            [+0.2631128296, +0.0083947774, +0.7284923930],
            [+0.2631128296, +0.7284923930, +0.0083947774],
            [+0.7284923930, +0.2631128296, +0.0083947774],
            [+0.7284923930, +0.0083947774, +0.2631128296],
            [+0.0083947774, +0.2631128296, +0.7284923930],
            [+0.0083947774, +0.7284923930, +0.2631128296],
        ]
        reference_weights = [
            +0.1443156077,
            +0.0950916343,
            +0.0950916343,
            +0.0950916343,
            +0.1032173705,
            +0.1032173705,
            +0.1032173705,
            +0.0324584976,
            +0.0324584976,
            +0.0324584976,
            +0.0272303142,
            +0.0272303142,
            +0.0272303142,
            +0.0272303142,
            +0.0272303142,
            +0.0272303142,
        ]
    else:
        raise ValueError("quadrature order not supported")
    return np.array(reference_points), np.array(reference_weights)
    # return reference_points, reference_weights


# def get_triangle_quadrature(
#     vertices: ndarray,
#     integration_order: int,
#     triangle_volume: float,
#     quadrature_type: QuadratureType = QuadratureType.GAUSS,
# ) -> (ndarray, ndarray):
#     """
#
#     :param vertices:
#     :param integration_order:
#     :param quadrature_type:
#     """
#     euclidean_dimension = vertices.shape[1]
#     q = Quadrature(integration_order, "TRIANGLE", quadrature_type=quadrature_type)
#     nq = Quadrature.get_number_of_quadrature_points(
#         integration_order, "TRIANGLE", quadrature_type=quadrature_type
#     )
#     if euclidean_dimension == 3:
#         p = get_triangle_reference_frame_transformation_matrix(vertices)
#         v0 = (p @ vertices.T).T
#         h = v0[0, 2]
#         quadrature_points = q.quadrature_table @ v0
#         quadrature_points[:, 2] = np.ones((nq,)) * h
#         p_inv = np.linalg.inv(p)
#         quadrature_points = (p_inv @ quadrature_points.T).T
#         # quadrature_weights = np.copy(q.quadrature_weights) * triangle_volume
#         quadrature_weights = (triangle_volume * np.eye(nq)) @ q.quadrature_weights
#     elif euclidean_dimension == 2:
#         quadrature_points = q.quadrature_table @ vertices
#         quadrature_weights = (triangle_volume * np.eye(nq)) @ q.quadrature_weights
#         # quadrature_weights = np.copy(q.quadrature_weights) * triangle_volume
#     else:
#         raise EnvironmentError("wrong")
#     return quadrature_points, quadrature_weights
