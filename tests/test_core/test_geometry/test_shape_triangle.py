from unittest import TestCase

from parameterized import parameterized
from scipy import integrate

from pythhon.parameters import *
from pythhon.geometry.shapes.segment import *
from pythhon.geometry.shapes.triangle import *
from pythhon.geometry.shape import Shape
from pythhon.fem.element.cell import Cell
from pythhon.fem.element.face import Face
from pythhon.fem.element.finite_element import FiniteElement


class Test(TestCase):
    def test_element_geometry(self):
        v0 = np.array([1.0, 1.3])
        v1 = np.array([2.1, 1.4])
        v2 = np.array([1.3, 1.6])
        # ---
        y0 = 1.0
        x0 = 1.2
        base = 2.1
        height = 0.6
        offset = 0.3
        theta = 0.2
        rotation_matrix = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
        v0 = np.array([x0, y0])
        v1 = np.array([x0 + base, y0])
        v2 = np.array([x0 + offset, y0 + height])
        triangle_vertices = np.zeros((2, 3))
        triangle_vertices[:, 0] = v0
        triangle_vertices[:, 1] = v1
        triangle_vertices[:, 2] = v2
        integration_order = 6
        rv = rotation_matrix @ triangle_vertices
        cell = Cell(ShapeType.TRIANGLE, rv, integration_order, QuadratureType.GAUSS)
        cell_centroid_check = np.array(
            [(rv[0, 0] + rv[0, 1] + rv[0, 2]) / 3.0, (rv[1, 0] + rv[1, 1] + rv[1, 2]) / 3.0]
        )
        faces = [
            Face(ShapeType.SEGMENT, rv[:,[0,1]], integration_order),
            Face(ShapeType.SEGMENT, rv[:,[1,2]], integration_order),
            Face(ShapeType.SEGMENT, rv[:,[2,0]], integration_order),
        ]

        # base = np.sqrt((v0[0] - v1[0]) ** 2 + (v0[1] - v1[1]) ** 2)
        cell_volume_check = (base * height) / 2.0
        cell_diameter_check = base
        np.testing.assert_allclose(cell.shape.centroid, cell_centroid_check, rtol=1e-10)
        np.testing.assert_allclose(cell.shape.volume, cell_volume_check, rtol=1e-10)
        np.testing.assert_allclose(cell.shape.diameter, cell_diameter_check, rtol=1e-10)
        self.assertTrue(True)
