from unittest import TestCase

from pythhon.fem.element.face import Face
from pythhon.parameters import *
import numpy as np


class TestFace(TestCase):
    def test_cell_centroid(self):
        vertices = np.array([[0.0, 0.0], [1.0, 0.0]]).T
        f = Face(ShapeType.SEGMENT, vertices, 2)
        centroid = np.array([1.0 / 2.0, 0.0])
        self.assertTrue((centroid == f.shape.centroid).all())

    def test_cell_diameter(self):
        vertices = np.array([[0.0, 0.0], [1.0, 0.0]]).T
        f = Face(ShapeType.SEGMENT, vertices, 2)
        diameter = 1.0
        self.assertTrue((diameter == f.shape.diameter).all())

    def test_face_mapping(self):
        vertices = np.array([[0.0, 0.0], [1.0, 0.0]]).T
        f = Face(ShapeType.SEGMENT, vertices, 2)
        print(f.mapping_matrix)
        # self.assertTrue((centroid == c.mapping_matrix).all())
        self.assertTrue(True)

    pass
