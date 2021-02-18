from unittest import TestCase

from parameterized import parameterized
from scipy import integrate

from pythhon.geometry.shapes.segment import *
from pythhon.math.linear_algebra import Vector, Matrix


class Test(TestCase):
    def test_get_segment_barycenter_cell(self):
        # vertices = np.array([[0.0], [1.0]])
        vertices = Matrix(1, 2, [[0.0], [1.0]])
        # vertices = geom.points(vertices)
        # barycenter = np.array([[1.0 / 2.0]])
        barycenter = Vector(1, [[1.0 / 2.0]])
        # barycenter_chek = get_segment_barycenter(vertices)
        test_passed = not (get_segment_barycenter(vertices) - barycenter).all()
        self.assertTrue(test_passed)

    def test_get_segment_barycenter_face(self):
        # vertices = np.array([[0.0, 0.0], [1.0, 0.0]])
        vertices = Matrix(2, 2, [[0.0, 0.0], [1.0, 0.0]])
        # vertices = geom.points(vertices)
        # barycenter = np.array([[1.0 / 2.0], [0.0]])
        barycenter = Vector(2, [[1.0 / 2.0], [0.0]])
        # barycenter_chek = get_segment_barycenter(vertices)
        test_passed = not (get_segment_barycenter(vertices) - barycenter).all()
        self.assertTrue(test_passed)

    def test_get_segment_volume_cell(self):
        vertices = np.array([[0.0], [1.0]])
        volume = 1.0
        # volume_check = get_segment_volume(vertices)
        # test_passed = not (get_segment_volume(vertices) - volume).all()
        test_passed = not get_segment_volume(vertices) - volume
        self.assertTrue(test_passed)

    def test_get_segment_volume_face(self):
        vertices = np.array([[0.0, 0.0], [1.0, 0.0]])
        volume = 1.0
        # volume_check = get_segment_volume(vertices)
        # test_passed = not (get_segment_volume(vertices) - volume).all()
        test_passed = not get_segment_volume(vertices) - volume
        self.assertTrue(test_passed)

    def test_get_segment_diameter_cell(self):
        vertices = np.array([[0.0], [1.2]])
        diameter = 1.2
        # volume_check = get_segment_volume(vertices)
        # test_passed = not (get_segment_volume(vertices) - volume).all()
        test_passed = not get_segment_diameter(vertices) - diameter
        self.assertTrue(test_passed)

    def test_get_segment_diameter_face(self):
        vertices = np.array([[0.0, 0.0], [1.2, 0.0]])
        diameter = 1.2
        # volume_check = get_segment_volume(vertices)
        # test_passed = not (get_segment_volume(vertices) - volume).all()
        test_passed = not get_segment_diameter(vertices) - diameter
        self.assertTrue(test_passed)

    # def test_get_segment_quadrature_cell(self):
    #     vertices = np.array([[0.0], [1.2]])
    #     # fun = (
    #     #     lambda y, x: 2.0 * x ** 3 * y ** 4
    #     #     + 1.0 * x ** 3
    #     #     - 3.0 * y ** 7
    #     #     + 1.0 * x
    #     #     - 6.0
    #     #     + 6.0 * y ** 7
    #     # )
    #     fun = (
    #         lambda x: 2.0 * x ** 8
    #         - 1.0 * x ** 3
    #         + 1.0 * x ** 7
    #         - 3.0 * x ** 2
    #         + 1.0 * x ** 1
    #         - 2.0
    #     )
    #     res = integrate.quad(fun, 0.0, 1.2)[0]
    #     quadrature_points, quadrature_weights = get_segment_quadrature(vertices, 8)
    #     res_num = np.sum(
    #         [fun(qp[0]) * qw for qp, qw in zip(quadrature_points, quadrature_weights)]
    #     )
    #     self.assertAlmostEqual(res, res_num)
    #
    # def test_get_segment_quadrature_face(self):
    #     vertices = np.array([[0.0, 0.0], [1.2, 0.0]])
    #     fun = (
    #         lambda x: 2.0 * x ** 8
    #         - 1.0 * x ** 3
    #         + 1.0 * x ** 7
    #         - 3.0 * x ** 2
    #         + 1.0 * x ** 1
    #         - 2.0
    #     )
    #     res = integrate.quad(fun, 0.0, 1.2)[0]
    #     quadrature_points, quadrature_weights = get_segment_quadrature(vertices, 8)
    #     res_num = np.sum(
    #         [fun(qp[0]) * qw for qp, qw in zip(quadrature_points, quadrature_weights)]
    #     )
    #     self.assertAlmostEqual(res, res_num)


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
