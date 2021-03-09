from typing import List

from pythhon.fem.element.elements.blocks import *
from pythhon.fem.element.finite_element import FiniteElement
from pythhon.pbbb.field import Field


def get_gradient_operators(field: Field, finite_element: FiniteElement, cell: Cell, faces: List[Face]):
    """

    Args:
        field:
        finite_element:
        cell:
        faces:

    Returns:

    """
    _d = field.euclidean_dimension
    _dx = field.field_dimension
    _cl = finite_element.cell_basis_l.dimension
    _ck = finite_element.cell_basis_k.dimension
    _fk = finite_element.face_basis_k.dimension
    _nf = len(faces)
    _es = _dx * (_cl + _nf * _fk)
    _gs = field.gradient_dimension
    nq = len(cell.quadrature_weights)
    gradient_operators = np.zeros((nq, _gs, _es), dtype=real)
    m_mas = np.zeros((_ck, _ck), dtype=real)
    for qc in range(len(cell.quadrature_weights)):
        x_q_c = cell.quadrature_points[:, qc]
        w_q_c = cell.quadrature_weights[qc]
        m_mas += get_cell_mass_matrix_in_cell(
            cell, finite_element.cell_basis_k, finite_element.cell_basis_k, x_q_c, w_q_c
        )
    m_mas_inv = np.linalg.inv(m_mas)
    for key, k in field.voigt_indices.items():
        i = key[0]
        j = key[1]
        local_grad_matric = np.zeros((_ck, _es), dtype=real)
        m_adv_j = np.zeros((_ck, _cl), dtype=real)
        for qc in range(len(cell.quadrature_weights)):
            x_q_c = cell.quadrature_points[:, qc]
            w_q_c = cell.quadrature_weights[qc]
            m_adv_j += get_cell_advection_matrix_in_cell(
                cell, finite_element.cell_basis_k, finite_element.cell_basis_l, j, x_q_c, w_q_c
            )
        _c0 = i * _cl
        _c1 = (i + 1) * _cl
        local_grad_matric[:, _c0:_c1] += m_adv_j
        for f, face in enumerate(faces):
            dist_in_face = (face.mapping_matrix @ (face.shape.centroid - cell.shape.centroid))[-1]
            # sign = dist[-1]
            if dist_in_face > 0:
                normal_vector_component = face.mapping_matrix[-1, j]
            else:
                normal_vector_component = -face.mapping_matrix[-1, j]
            m_mas_f = np.zeros((_ck, _cl), dtype=real)
            m_hyb_f = np.zeros((_ck, _fk), dtype=real)
            for qf in range(len(face.quadrature_weights)):
                x_q_f = face.quadrature_points[:, qf]
                w_q_f = face.quadrature_weights[qf]
                m_mas_f += get_cell_mass_matrix_in_face(
                    cell, finite_element.cell_basis_k, finite_element.cell_basis_l, x_q_f, w_q_f
                )
                m_hyb_f += get_hybrid_mass_matrix_in_face(
                    cell, face, finite_element.cell_basis_k, finite_element.face_basis_k, x_q_f, w_q_f
                )
            _c0 = i * _cl
            _c1 = (i + 1) * _cl
            local_grad_matric[:, _c0:_c1] -= m_mas_f * normal_vector_component
            _c0 = _dx * _cl + f * _dx * _fk + i * _fk
            _c1 = _dx * _cl + f * _dx * _fk + (i + 1) * _fk
            local_grad_matric[:, _c0:_c1] += m_hyb_f * normal_vector_component
        local_grad_matric = m_mas_inv @ local_grad_matric
        for qc in range(len(cell.quadrature_weights)):
            x_q_c = cell.quadrature_points[:, qc]
            v = finite_element.cell_basis_k.evaluate_function(x_q_c, cell.shape.centroid, cell.shape.diameter)
            gradient_component = field.voigt_coefficients[key] * v @ local_grad_matric
            gradient_operators[qc, k] += gradient_component
    return gradient_operators