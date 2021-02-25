from unittest import TestCase

from parameterized import parameterized
from scipy import integrate

from pythhon.parameters import *
from pythhon.geometry.shapes.segment import *
from pythhon.geometry.shape import Shape

class TestSequence(TestCase):
    @parameterized.expand(
        [
            [(lambda x: 2.0 * x ** 1 + 3.0), 1],
            [(lambda x: 1.0 * x ** 2 + 6.0 * x ** 1 - 2.0), 2],
            [(lambda x: 4.0 * x ** 3 - 4.0 * x ** 2 + 5.0 * x ** 1 - 1.0), 3],
            [(lambda x: 7.0 * x ** 4 - 4.0 * x ** 3 + 5.0 * x ** 2 - 1.0 * x ** 1 + 9.0), 4,],
            [(lambda x: 2.0 * x ** 5 - 1.0 * x ** 4 + 1.0 * x ** 3 - 3.0 * x ** 2 + 1.0 * x ** 1 - 6.0), 5,],
            [(lambda x: 2.0 * x ** 6 - 1.0 * x ** 4 + 1.0 * x ** 3 - 3.0 * x ** 2 + 1.0 * x ** 1 - 9.0), 6,],
            [(lambda x: 2.0 * x ** 7 - 1.0 * x ** 5 + 1.0 * x ** 6 - 3.0 * x ** 2 + 1.0 * x ** 1 - 4.0), 7,],
            [(lambda x: 2.0 * x ** 8 - 1.0 * x ** 3 + 1.0 * x ** 7 - 3.0 * x ** 2 + 1.0 * x ** 1 - 2.0), 8,],
        ]
    )
    def test_get_segment_quadrature_cell(self, test_polynomial, integration_order):
        vertices = np.array([[0.0], [1.2]])
        res = integrate.quad(test_polynomial, 0.0, 1.2)[0]
        # res = integrate.dblquad(
        #     test_polynomial, 0.0, 1.2, lambda x: 0.0, lambda x: 1 - (x / 1.2)
        # )[0]
        quadrature_points, quadrature_weights = get_segment_quadrature(vertices, integration_order)
        res_num = np.sum([test_polynomial(qp[0]) * qw for qp, qw in zip(quadrature_points, quadrature_weights)])
        self.assertAlmostEqual(res, res_num)


class TestSequence(TestCase):
    @parameterized.expand(
        [
            [(lambda x: 2.0 * x ** 1 + 3.0), 1],
            [(lambda x: 1.0 * x ** 2 + 6.0 * x ** 1 - 2.0), 2],
            [(lambda x: 4.0 * x ** 3 - 4.0 * x ** 2 + 5.0 * x ** 1 - 1.0), 3],
            [(lambda x: 7.0 * x ** 4 - 4.0 * x ** 3 + 5.0 * x ** 2 - 1.0 * x ** 1 + 9.0), 4,],
            [(lambda x: 2.0 * x ** 5 - 1.0 * x ** 4 + 1.0 * x ** 3 - 3.0 * x ** 2 + 1.0 * x ** 1 - 6.0), 5,],
            [(lambda x: 2.0 * x ** 6 - 1.0 * x ** 4 + 1.0 * x ** 3 - 3.0 * x ** 2 + 1.0 * x ** 1 - 9.0), 6,],
            [(lambda x: 2.0 * x ** 7 - 1.0 * x ** 5 + 1.0 * x ** 6 - 3.0 * x ** 2 + 1.0 * x ** 1 - 4.0), 7,],
            [(lambda x: 2.0 * x ** 8 - 1.0 * x ** 3 + 1.0 * x ** 7 - 3.0 * x ** 2 + 1.0 * x ** 1 - 2.0), 8,],
        ]
    )
    def test_get_segment_quadrature_face(self, test_polynomial, integration_order):
        vertices = np.array([[0.0, 0.0], [1.2, 0.0]])
        # res = integrate.dblquad(
        #     test_polynomial, 0.0, 1.2, lambda x: 0.0, lambda x: 1 - (x / 1.2)
        # )[0]
        res = integrate.quad(test_polynomial, 0.0, 1.2)[0]
        quadrature_points, quadrature_weights = get_segment_quadrature(vertices, integration_order)
        res_num = np.sum([test_polynomial(qp[0]) * qw for qp, qw in zip(quadrature_points, quadrature_weights)])
        self.assertAlmostEqual(res, res_num)