from typing import List

from numpy import ndarray
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

from pythhon.fem.element.finite_element import FiniteElement
from pythhon.mesh.mesh import Mesh
from pythhon.parameters import *
from pythhon.pbbb.boundary_condition import BoundaryCondition
from pythhon.pbbb.field import Field


def get_total_system_size(
    field: Field, finite_element: FiniteElement, mesh: Mesh, boundary_conditions: List[BoundaryCondition]
) -> (int, int):
    """

    Args:
        finite_element:
        mesh:
        boundary_conditions:

    Returns:

    """
    constrained_faces = 0
    constrained_constants = 0
    for key, val in mesh.faces_boundaries_connectivity.items():
        for bc in boundary_conditions:
            if key == bc.boundary_name and bc.boundary_type == BoundaryType.DISPLACEMENT:
                constrained_faces += len(val)
            elif key == bc.boundary_name and bc.boundary_type == BoundaryType.SLIDE:
                constrained_constants += len(val)
    constrained_faces_matrix_size = constrained_faces * finite_element.face_basis_k.dimension
    constrained_const_matrix_size = constrained_constants
    lagrange_system_size = constrained_faces_matrix_size + constrained_const_matrix_size
    system_size = mesh.number_of_faces_in_mesh * finite_element.face_basis_k.dimension * field.field_dimension
    constrained_system_size = system_size + lagrange_system_size
    return constrained_system_size, system_size


def solve_system(tangent_matrix: ndarray, residual: ndarray) -> ndarray:
    """

    Args:
        tangent_matrix:
        residual:

    Returns:

    """
    sparse_global_matrix = csr_matrix(-tangent_matrix)
    correction = spsolve(sparse_global_matrix, residual)
    return correction
