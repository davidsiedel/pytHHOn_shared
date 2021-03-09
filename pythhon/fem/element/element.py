from numpy import ndarray
from typing import List, Dict, Tuple, Any
from pythhon.parameters import *
from pythhon.pbbb.boundary_condition import BoundaryCondition
from pythhon.pbbb.load import Load
from pythhon.fem.element.cell import Cell
from pythhon.fem.element.face import Face
from pythhon.fem.element.finite_element import FiniteElement
from pythhon.pbbb.field import Field

import pythhon.fem.element.elements.hdgs as hdgs
import pythhon.fem.element.elements.hdg as hdg
import pythhon.fem.element.elements.stabilization as stab


class Element:
    cell: Cell
    faces: List[Face]
    gradients_operators: ndarray
    stabilization_operator: ndarray
    element_size: int
    faces_indices: List[int]
    faces_lagrange_system_row_positions: List[List[Tuple[Any, Any]]]
    faces_lagrange_local_positions: List[List[Tuple[Any, Any]]]
    m_cell_cell_inv: ndarray
    m_cell_faces: ndarray
    v_cell: ndarray
    cell_unknown_vector: ndarray

    def __init__(
        self,
        field: Field,
        finite_element: FiniteElement,
        cell: Cell,
        faces: List[Face],
        faces_indices: List[int],
        faces_lagrange_system_row_positions: List[List[Tuple[Any, Any]]],
        faces_lagrange_system_col_positions: List[List[Tuple[Any, Any]]],
        faces_lagrange_local_positions: List[List[Tuple[Any, Any]]],
    ):
        """

        Args:
            field:
            finite_element:
            cell:
            faces:
            faces_indices:
            faces_lagrange_positions:
            faces_lag_local_positions:
        """
        self.cell = cell
        self.faces = faces
        self.faces_indices = faces_indices
        self.faces_lagrange_system_row_positions = faces_lagrange_system_row_positions
        self.faces_lagrange_system_col_positions = faces_lagrange_system_col_positions
        self.faces_lagrange_local_positions = faces_lagrange_local_positions
        _dx = field.field_dimension
        _cl = finite_element.cell_basis_l.dimension
        _fk = finite_element.face_basis_k.dimension
        _nf = len(self.faces)
        self.element_size = _dx * (_cl + _nf * _fk)
        self.m_cell_cell_inv = np.zeros((_cl * _dx, _cl * _dx), dtype=real)
        self.m_cell_faces = np.zeros((_cl * _dx, _nf * _fk * _dx), dtype=real)
        self.v_cell = np.zeros((_cl * _dx,), dtype=real)
        if finite_element.element_type in [ElementType.HDG_LOW, ElementType.HDG_EQUAL, ElementType.HDG_HIGH]:
            if field.derivation_type == DerivationType.SYMMETRIC:
                self.gradients_operators = hdgs.get_gradient_operators(field, finite_element, cell, faces)
                # self.stabilization_operator = hdgs.get_stabilization_operator(field, finite_element, cell, faces)
                self.stabilization_operator = stab.get_stab_test(field, finite_element, cell, faces)
                # self.stabilization_operator = stab.get_stabilization_operator(field, finite_element, cell, faces)
            elif field.derivation_type == DerivationType.FULL:
                self.gradients_operators = hdg.get_gradient_operators(field, finite_element, cell, faces)
                # self.stabilization_operator = hdg.get_stabilization_operator(field, finite_element, cell, faces)
                self.stabilization_operator = stab.get_stab_test(field, finite_element, cell, faces)
                # self.stabilization_operator = stab.get_stabilization_operator(field, finite_element, cell, faces)
            else:
                raise KeyError("NO")
        else:
            raise KeyError("NO")
        self.cell_unknown_vector = np.zeros((_dx * _cl,), dtype=real)
        return

    def get_element_unknown_vector(
        self, field: Field, finite_element: FiniteElement, faces_unknown_vector: ndarray
    ) -> ndarray:
        """

        Args:
            field:
            finite_element:
            faces_unknown_vector:

        Returns:

        """
        _dx = field.field_dimension
        _fk = finite_element.face_basis_k.dimension
        _cl = finite_element.cell_basis_l.dimension
        element_unknown_vector = np.zeros((self.element_size,), dtype=real)
        for _i_local, _i_global in enumerate(self.faces_indices):
            _c0_fg = _i_global * (_fk * _dx)
            _c1_fg = (_i_global + 1) * (_fk * _dx)
            _c0_fl = (_dx * _cl) + _i_local * (_fk * _dx)
            _c1_fl = (_dx * _cl) + (_i_local + 1) * (_fk * _dx)
            element_unknown_vector[_c0_fl:_c1_fl] += faces_unknown_vector[_c0_fg:_c1_fg]
        _c0_c = _cl * _dx
        element_unknown_vector[:_c0_c] += self.cell_unknown_vector
        return element_unknown_vector

    def get_transformation_gradient(
        self, field: Field, finite_element: FiniteElement, faces_unknown_vector: ndarray, _qc: int
    ):
        """

        Args:
            field:
            finite_element:
            faces_unknown_vector:
            _qc:

        Returns:

        """
        element_unknown_vector = self.get_element_unknown_vector(field, finite_element, faces_unknown_vector)
        if field.strain_type == StrainType.DISPLACEMENT_TRANSFORMATION_GRADIENT:
            ones_vect = np.zeros((field.gradient_dimension,), dtype=real)
            for indx in range(3):
                ones_vect[indx] += 1.0
            # transformation_gradient = ones_vect + self.gradients_operators[_qc] @ element_unknown_vector
            transformation_gradient_0 = self.gradients_operators[_qc] @ element_unknown_vector
            # transformation_gradient_0[1:] = np.zeros((field.gradient_dimension-1,))
            transformation_gradient = ones_vect + transformation_gradient_0
        elif field.strain_type == StrainType.DISPLACEMENT_SYMMETRIC_GRADIENT:
            # transformation_gradient = self.gradients_operators[_qc] @ element_unknown_vector
            transformation_gradient_0 = self.gradients_operators[_qc] @ element_unknown_vector
            # transformation_gradient_0[-1] = 0.0
            # transformation_gradient_0[1] = 0.0
            transformation_gradient = transformation_gradient_0
        else:
            raise KeyError("No such strain type")
        return transformation_gradient

    def get_cell_field_value(
        self, field: Field, finite_element: FiniteElement, faces_unknown_vector: ndarray, point: ndarray, direction: int
    ) -> float:
        """

        Args:
            field:
            finite_element:
            faces_unknown_vector:
            point:
            direction:

        Returns:

        """
        element_unknown_vector = self.get_element_unknown_vector(field, finite_element, faces_unknown_vector)
        _cl = finite_element.cell_basis_l.dimension
        _c0 = direction * _cl
        _c1 = (direction + 1) * _cl
        v = finite_element.cell_basis_l.evaluate_function(point, self.cell.shape.centroid, self.cell.shape.diameter)
        field_unknown_value = v @ element_unknown_vector[_c0:_c1]
        return field_unknown_value