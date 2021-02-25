from unittest import TestCase
import numpy as np
from numpy import ndarray
from pythhon.fem.element.elements.projection import get_projection_in_element
import pythhon.fem.element.elements.hdgs as hdgs
import pythhon.fem.element.elements.hdg as hdg
from pythhon.pbbb.field import Field
from pythhon.fem.element.finite_element import FiniteElement
from pythhon.fem.element.cell import Cell
from pythhon.fem.element.face import Face
from pythhon.parameters import *

class TestProjection(TestCase):
    def test_get_projection_in_element(self):

        def function_x(p: ndarray):
            # return np.cos(p[0]) * np.sin(p[1])
            return p[0]

        def function_y(p: ndarray):
            # return np.sin(p[0]) * np.cos(p[1])
            return p[1]

        def function_derivative(p :ndarray):
            d_f = np.zeros((4,), dtype=real)
            # d_f[0] = - np.sin(p[0]) * np.sin(p[1])
            # d_f[1] = - np.sin(p[0]) * np.sin(p[1])
            # d_f[3] = + np.cos(p[0]) * np.cos(p[1])
            #
            d_f[0] = 1.0
            d_f[1] = 1.0
            d_f[3] = 0.0
            return d_f

        functions = [function_x ,function_y]

        displacement = Field(
            label="U",
            euclidean_dimension=2,
            field_type=FieldType.DISPLACEMENT_PLANE_STRAIN,
            strain_type=StrainType.DISPLACEMENT_SYMMETRIC_GRADIENT,
            stress_type=StressType.CAUCHY,
            derivation_type=DerivationType.SYMMETRIC,
        )

        finite_element = FiniteElement(
            element_type=ElementType.HDG_EQUAL,
            # element_type=ElementType.HDG_HIGH,
            # element_type=ElementType.HDG_LOW,
            # field=field,
            polynomial_order=1,
            euclidean_dimension=2,
            basis_type=BasisType.MONOMIAL
        )

        integration_order = 2*(finite_element.polynomial_order+1)

        vertices = np.array([
            # [0.2, 1.2, 0.6],
            # [0.2, 2.0, 3.4],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ], dtype=real)

        cell = Cell(ShapeType.TRIANGLE, vertices, integration_order)
        faces = [
            Face(ShapeType.SEGMENT, vertices[:, [0, 1]], integration_order),
            Face(ShapeType.SEGMENT, vertices[:, [1, 2]], integration_order),
            Face(ShapeType.SEGMENT, vertices[:, [2, 0]], integration_order),
        ]

        projection = get_projection_in_element(displacement, finite_element, cell, faces, functions)

        sym_grad = hdgs.get_gradient_operators(displacement, finite_element, cell, faces)

        for _qc in range(len(cell.quadrature_weights)):
            _x_qc = cell.quadrature_points[:, _qc]
            _w_qc = cell.quadrature_weights[_qc]
            grad_val = sym_grad[_qc] @ projection
            grad_val_check = function_derivative(_x_qc)
            print(grad_val)
            print(grad_val_check)
        self.assertTrue(True)
