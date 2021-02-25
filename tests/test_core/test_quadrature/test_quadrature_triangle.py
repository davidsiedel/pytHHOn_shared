from unittest import TestCase

from parameterized import parameterized
from scipy import integrate

from pythhon.parameters import *
from pythhon.geometry.shapes.segment import *
from pythhon.geometry.shape import Shape

class TestSequence(TestCase):
    @parameterized.expand(
        [
            [(lambda y, x: 2.0 * x ** 1 + 3.0 * y ** 1 + 7.0), 1],
            [(lambda y, x: 1.0 * x ** 2 + 6.0 * x + 3.0 * y ** 2 - 2.0 * y + 2.0 * x * y - 2.0), 2,],
            [(lambda y, x: 4.0 * x ** 3 - 4.0 * x ** 2 + 5.0 * x - 1.0 + 4.0 * y ** 3), 3,],
            [(lambda y, x: 7.0 * x ** 4 - 4.0 * x ** 3 + 5.0 * x ** 2 - 1.0 * x + 9.0 + 5.0 * y ** 4), 3,],
            [
                (lambda y, x: 2.0 * x ** 5 - 1.0 * x ** 4 + 1.0 * x ** 3 - 3.0 * x ** 2 + 1.0 * x - 6.0 + 6.0 * y ** 5),
                4,
            ],
            [(lambda y, x: 2.0 * x ** 3 * y ** 3 + 1.0 * x ** 3 - 3.0 * x ** 2 + 1.0 * x - 6.0 + 6.0 * y ** 6), 5,],
            [(lambda y, x: 2.0 * x ** 3 * y ** 4 + 1.0 * x ** 3 - 3.0 * y ** 7 + 1.0 * x - 6.0 + 6.0 * y ** 7), 6,],
            [(lambda y, x: 2.0 * x ** 4 * y ** 4 + 1.0 * x ** 3 - 3.0 * y ** 7 + 1.0 * x - 6.0 + 6.0 * y ** 7), 7,],
            [
                (
                    lambda y, x: 2.0 * (x ** 4) * (y ** 4)
                    + 1.0 * x ** 3
                    - 3.0 * y ** 7
                    + 1.0 * x
                    - 6.0
                    + 6.0 * y ** 8
                    - 123.0 * x ** 8
                    + 530.0 * (x ** 3) * (y ** 4)
                ),
                8,
            ],
        ]
    )
    def test_get_triangle_quadrature_cell(self, test_polynomial, integration_order):
        vertices = np.array([[0.0, 0.0], [1.2, 0.0], [0.0, 1.0]])
        res = integrate.dblquad(test_polynomial, 0.0, 1.2, lambda x: 0.0, lambda x: 1 - (x / 1.2))[0]
        quadrature_points, quadrature_weights = get_triangle_quadrature(vertices, integration_order)
        res_num = np.sum([test_polynomial(qp[1], qp[0]) * qw for qp, qw in zip(quadrature_points, quadrature_weights)])
        self.assertAlmostEqual(res, res_num)


class TestSequence(TestCase):
    @parameterized.expand(
        [
            [(lambda y, x: 2.0 * x ** 1 + 3.0 * y ** 1 + 7.0), 1],
            [(lambda y, x: 1.0 * x ** 2 + 6.0 * x + 3.0 * y ** 2 - 2.0 * y + 2.0 * x * y - 2.0), 2,],
            [(lambda y, x: 4.0 * x ** 3 - 4.0 * x ** 2 + 5.0 * x - 1.0 + 4.0 * y ** 3), 3,],
            [(lambda y, x: 7.0 * x ** 4 - 4.0 * x ** 3 + 5.0 * x ** 2 - 1.0 * x + 9.0 + 5.0 * y ** 4), 3,],
            [
                (lambda y, x: 2.0 * x ** 5 - 1.0 * x ** 4 + 1.0 * x ** 3 - 3.0 * x ** 2 + 1.0 * x - 6.0 + 6.0 * y ** 5),
                4,
            ],
            [(lambda y, x: 2.0 * x ** 3 * y ** 3 + 1.0 * x ** 3 - 3.0 * x ** 2 + 1.0 * x - 6.0 + 6.0 * y ** 6), 5,],
            [(lambda y, x: 2.0 * x ** 3 * y ** 4 + 1.0 * x ** 3 - 3.0 * y ** 7 + 1.0 * x - 6.0 + 6.0 * y ** 7), 6,],
            [(lambda y, x: 2.0 * x ** 4 * y ** 4 + 1.0 * x ** 3 - 3.0 * y ** 7 + 1.0 * x - 6.0 + 6.0 * y ** 7), 7,],
            [
                (
                    lambda y, x: 2.0 * (x ** 4) * (y ** 4)
                    + 1.0 * x ** 3
                    - 3.0 * y ** 7
                    + 1.0 * x
                    - 6.0
                    + 6.0 * y ** 8
                    - 123.0 * x ** 8
                    + 530.0 * (x ** 3) * (y ** 4)
                ),
                8,
            ],
        ]
    )
    def test_get_triangle_quadrature_face(self, test_polynomial, integration_order):
        vertices = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.0, 1.0, 0.0]])
        res = integrate.dblquad(test_polynomial, 0.0, 1.2, lambda x: 0.0, lambda x: 1 - (x / 1.2))[0]
        quadrature_points, quadrature_weights = get_triangle_quadrature(vertices, integration_order)
        res_num = np.sum([test_polynomial(qp[1], qp[0]) * qw for qp, qw in zip(quadrature_points, quadrature_weights)])
        self.assertAlmostEqual(res, res_num)