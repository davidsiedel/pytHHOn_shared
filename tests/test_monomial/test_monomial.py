from unittest import TestCase

from parameterized import parameterized

from pythhon.fem.basis.bases.monomial import *


class TestSequence(TestCase):
    @parameterized.expand(
        [
            [np.array([[0]]), 0, 0],
            [np.array([[0], [1], [2]]), 2, 1],
            [np.array([[0, 0], [0, 1], [1, 0], [0, 2], [1, 1], [2, 0]]), 2, 2],
            [np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0]]), 1, 3],
        ]
    )
    def test_get_exponents(self, exponents, polynomial_order, euclidean_dimension):
        passed = not (get_exponents(polynomial_order, euclidean_dimension) - exponents).all()
        self.assertTrue(passed)


class TestSequence(TestCase):
    @parameterized.expand(
        [
            [np.array([[0]]), 0, 0],
            [np.array([[0], [1], [2]]), 2, 1],
            [np.array([[0, 0], [0, 1], [1, 0], [0, 2], [1, 1], [2, 0]]), 2, 2],
            [np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0]]), 1, 3],
        ]
    )
    def test_get_exponents(self, exponents, polynomial_order, euclidean_dimension):
        passed = not (get_exponents(polynomial_order, euclidean_dimension) - exponents).all()
        self.assertTrue(passed)
