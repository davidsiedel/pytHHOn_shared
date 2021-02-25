from unittest import TestCase

from pythhon.fem.element.cell import Cell
from pythhon.parameters import *
import numpy as np


class TestCell(TestCase):
    # cell_shape_type: ShapeType,
    # cell_vertices: ndarray,
    # integration_order: int,
    # quadrature_type: QuadratureType = QuadratureType.GAUSS,

    def test_cell_centroid(self):
        vertices = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]).T
        c = Cell(ShapeType.TRIANGLE, vertices, 2)
        centroid = np.array([1.0 / 3.0, 1.0 / 3.0])
        self.assertTrue((centroid == c.shape.centroid).all())

    def test_cell_diameter(self):
        vertices = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]).T
        c = Cell(ShapeType.TRIANGLE, vertices, 2)
        diameter = np.sqrt(2.0)
        print(c.quadrature_points)
        print(c.quadrature_weights)
        self.assertTrue((diameter == c.shape.diameter).all())

    pass
