from pythhon.parameters import *
from pythhon.fem.basis.basis import Basis
from pythhon.pbbb.field import Field
from typing import Dict, Tuple


class FiniteElement:
    element_type: ElementType
    polynomial_order: int
    basis_type: BasisType
    cell_basis_k: Basis
    cell_basis_l: Basis
    face_basis_k: Basis
    cell_basis_r: Basis

    def __init__(
        self,
        element_type: ElementType,
        polynomial_order: int,
        euclidean_dimension: int,
        basis_type: BasisType = BasisType.MONOMIAL,
    ):
        """

        Args:
            element_type:
            polynomial_order:
            euclidean_dimension:
            basis_type:
        """
        self.element_type = element_type
        self.polynomial_order = polynomial_order
        self.basis_type = basis_type
        if element_type == ElementType.HDG_EQUAL:
            self.cell_basis_k = Basis(
                polynomial_order, euclidean_dimension, basis_type=basis_type
            )
            self.cell_basis_l = Basis(
                polynomial_order, euclidean_dimension, basis_type=basis_type
            )
            self.face_basis_k = Basis(
                polynomial_order, euclidean_dimension - 1, basis_type=basis_type
            )
            self.cell_basis_r = Basis(
                polynomial_order + 1, euclidean_dimension, basis_type=basis_type
            )
        elif element_type == ElementType.HDG_HIGH:
            self.cell_basis_k = Basis(
                polynomial_order, euclidean_dimension, basis_type=basis_type
            )
            self.cell_basis_l = Basis(
                polynomial_order + 1, euclidean_dimension, basis_type=basis_type
            )
            self.face_basis_k = Basis(
                polynomial_order, euclidean_dimension - 1, basis_type=basis_type
            )
            self.cell_basis_r = Basis(
                polynomial_order + 1, euclidean_dimension, basis_type=basis_type
            )
        elif element_type == ElementType.HDG_LOW:
            self.cell_basis_k = Basis(
                polynomial_order, euclidean_dimension, basis_type=basis_type
            )
            self.cell_basis_l = Basis(
                polynomial_order - 1, euclidean_dimension, basis_type=basis_type
            )
            self.face_basis_k = Basis(
                polynomial_order, euclidean_dimension - 1, basis_type=basis_type
            )
            self.cell_basis_r = Basis(
                polynomial_order + 1, euclidean_dimension, basis_type=basis_type
            )
        else:
            raise KeyError("NO")
        return
