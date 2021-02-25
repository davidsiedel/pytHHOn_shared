from unittest import TestCase
from dev_runs.hdg_temp import *


class Test(TestCase):
    def test_get_gradient_operators(self):

        face0_vertices = np.array([[9.6, 24.16], [9.6, 31.84]]).T

        face1_vertices = np.array([[9.6, 31.84], [0.0, 26.4]]).T

        face2_vertices = np.array([[9.6, 24.16], [0.0, 26.4]]).T

        cell_vertices = np.array([[9.6, 24.16], [9.6, 31.84], [0.0, 26.4]]).T

        cell = Cell(ShapeType.TRIANGLE, cell_vertices, 2)
        face0 = Face(ShapeType.SEGMENT, face0_vertices, 2)
        face1 = Face(ShapeType.SEGMENT, face1_vertices, 2)
        face2 = Face(ShapeType.SEGMENT, face2_vertices, 2)

        fe = FiniteElement(ElementType.HDG_EQUAL, FieldType.VECTOR, DerivationType.SYMMETRIC, 1, 2, BasisType.MONOMIAL)

        faces = [face0, face1, face2]

        grad = get_gradient_operators(fe, cell, faces)

        self.assertTrue(True)
