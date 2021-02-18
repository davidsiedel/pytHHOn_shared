from typing import List
from numpy import ndarray
import matplotlib.pyplot as plt

from pythhon.mesh.mesh import Mesh
from pythhon.fem.element.cell import Cell
from pythhon.fem.element.face import Face
from pythhon.fem.element.element import Element
from pythhon.fem.element.finite_element import FiniteElement
from pythhon.pbbb.boundary_condition import BoundaryCondition
from pythhon.pbbb.load import Load
from pythhon.pbbb.material import Material
from pythhon.pbbb.field import Field
import pythhon.fem.element.elements.blocks as blocks
from pythhon.parameters import *

import pythhon.pbbb.system as systm
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix

from mgis import behaviour as mgis_bv


def get_integration_order(element_type: ElementType, polynomial_order: intg) -> intg:
    """

    Args:
        element_type:
        polynomial_order:

    Returns:

    """
    if element_type in [ElementType.HHO_LOW, ElementType.HDG_LOW]:
        integration_order = 2 * (polynomial_order + 1)
    elif element_type in [ElementType.HHO_EQUAL, ElementType.HDG_EQUAL]:
        integration_order = 2 * (polynomial_order + 1)
    elif element_type in [ElementType.HHO_HIGH, ElementType.HDG_HIGH]:
        integration_order = 2 * (polynomial_order + 1)
    else:
        raise KeyError("NO")
    return integration_order


class Problem:
    mesh: Mesh
    finite_element: FiniteElement
    elements: List[Element]

    def __init__(
        self,
        mesh_file_path: str,
        field: Field,
        polynomial_order: int,
        time_steps: ndarray,
        iterations: int,
        finite_element: FiniteElement,
        boundary_conditions: List[BoundaryCondition],
        loads: List[Load] = None,
        quadrature_type: QuadratureType = QuadratureType.GAUSS,
        tolerance: float = 1.0e-6,
    ):
        """

        Args:
            mesh_file_path:
            field:
            polynomial_order:
            time_steps:
            iterations:
            finite_element:
            boundary_conditions:
            loads:
            quadrature_type:
            tolerance:
        """
        integration_order = get_integration_order(finite_element.element_type, polynomial_order)
        self.mesh = Mesh(mesh_file_path=mesh_file_path, integration_order=integration_order)
        self.field = field
        self.finite_element = finite_element
        self.__check_loads(loads)
        self.__check_boundary_conditions(boundary_conditions)
        self.boundary_conditions = boundary_conditions
        self.loads = loads
        self.time_steps = time_steps
        self.number_of_iterations = iterations
        self.tolerance = tolerance
        # ------ build elements
        self.elements = self.get_elements(integration_order, quadrature_type=quadrature_type)
        return

    def clean_res_dir(self):
        """

        Returns:

        """
        res_folder = os.path.join(get_project_path(), "res")
        print(res_folder)
        for filename in os.listdir(res_folder):
            file_path = os.path.join(res_folder, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print("Failed to delete %s. Reason: %s" % (file_path, e))
        return

    def create_vertex_res_files(self, suffix: str):
        """

        Args:
            suffix:

        Returns:

        """
        with open(os.path.join(get_project_path(), "res/res_vtx_{}.csv".format(suffix)), "w") as res_vtx_file:
            for x_dir in range(self.field.euclidean_dimension):
                res_vtx_file.write("X_{},".format(x_dir))
            for u_dir in range(self.field.field_dimension):
                res_vtx_file.write("{}_{},".format(self.field.label, u_dir))
            res_vtx_file.write("\n")
        return

    def create_quadrature_points_res_files(self, suffix: str):
        """

        Args:
            suffix:
        """
        with open(os.path.join(get_project_path(), "res/res_qdp_{}.csv".format(suffix)), "w") as res_qdp_file:
            for x_dir in range(self.field.euclidean_dimension):
                res_qdp_file.write("XQ_{},".format(x_dir))
            for u_dir in range(self.field.field_dimension):
                res_qdp_file.write("{}_{},".format(self.field.label, u_dir))
            for strain_component in range(self.field.gradient_dimension):
                res_qdp_file.write("STRAIN_{},".format(strain_component))
            for stress_component in range(self.field.gradient_dimension):
                res_qdp_file.write("STRESS_{},".format(stress_component))
            res_qdp_file.write("\n")

    def write_vertex_res_files(self, suffix: str, unknown_increment: ndarray):
        """

        Args:
            suffix:
            unknown_increment:

        Returns:

        """
        with open(os.path.join(get_project_path(), "res/res_vtx_{}.csv".format(suffix)), "a") as res_vtx_file:
            for vertex_count in range(self.mesh.number_of_vertices_in_mesh):
                vertex = self.mesh.vertices[:, vertex_count]
                for x_dir in range(self.field.euclidean_dimension):
                    res_vtx_file.write("{},".format(vertex[x_dir]))
                vertex_field_value = np.zeros((self.field.field_dimension,), dtype=real)
                for c, cell_vertices_connectivity in enumerate(self.mesh.cells_vertices_connectivity):
                    if vertex_count in cell_vertices_connectivity:
                        for u_dir in range(self.field.field_dimension):
                            vertex_field_value[u_dir] += self.elements[c].get_cell_field_increment_value(
                                point=vertex,
                                direction=u_dir,
                                field=self.field,
                                finite_element=self.finite_element,
                                global_vector=unknown_increment,
                            )
                vertex_field_value = vertex_field_value / self.mesh.vertices_weights_cell[vertex_count]
                for u_dir in range(self.field.field_dimension):
                    res_vtx_file.write("{},".format(vertex_field_value[u_dir]))
                res_vtx_file.write("\n")
        return

    def write_quadrature_points_res_files(self, suffix: str, unknown_increment: ndarray, material: Material):
        """

        Args:
            suffix:
            unknown_increment:
            material:

        Returns:

        """
        with open(os.path.join(get_project_path(), "res/res_qdp_{}.csv".format(suffix)), "a") as res_qdp_file:
            qp = 0
            for element in self.elements:
                for qc in range(len(element.cell.quadrature_weights)):
                    x_q_c = element.cell.quadrature_points[:, qc]
                    for x_dir in range(self.field.euclidean_dimension):
                        res_qdp_file.write("{},".format(x_q_c[x_dir]))
                    for u_dir in range(self.field.field_dimension):
                        quad_point_field_value = element.get_cell_field_increment_value(
                            point=x_q_c,
                            direction=u_dir,
                            field=self.field,
                            finite_element=self.finite_element,
                            global_vector=unknown_increment,
                        )
                        res_qdp_file.write("{},".format(quad_point_field_value))
                    for g_dir in range(self.field.gradient_dimension):
                        strain_component = material.mat_data.s1.gradients[qp, g_dir]
                        res_qdp_file.write("{},".format(strain_component))
                    for g_dir in range(self.field.gradient_dimension):
                        if self.field.strain_type == StrainType.DISPLACEMENT_TRANSFORMATION_GRADIENT:
                            F = np.zeros((3, 3))
                            F[0, 0] = material.mat_data.s1.gradients[qp, 0]
                            F[1, 1] = material.mat_data.s1.gradients[qp, 1]
                            F[2, 2] = material.mat_data.s1.gradients[qp, 2]
                            F[0, 1] = material.mat_data.s1.gradients[qp, 3]
                            F[1, 0] = material.mat_data.s1.gradients[qp, 4]
                            PK = np.zeros((3, 3))
                            PK[0, 0] = material.mat_data.s1.thermodynamic_forces[qp, 0]
                            PK[1, 1] = material.mat_data.s1.thermodynamic_forces[qp, 1]
                            PK[2, 2] = material.mat_data.s1.thermodynamic_forces[qp, 2]
                            PK[0, 1] = material.mat_data.s1.thermodynamic_forces[qp, 3]
                            PK[1, 0] = material.mat_data.s1.thermodynamic_forces[qp, 4]
                            J = np.linalg.det(F)
                            # F_T_inv = np.linalg.inv(F.T)
                            sig = (1.0 / J) * PK @ F.T
                            sig_vect = np.zeros((5,))
                            sig_vect[0] = sig[0, 0]
                            sig_vect[1] = sig[1, 1]
                            sig_vect[2] = sig[2, 2]
                            sig_vect[3] = sig[0, 1]
                            sig_vect[4] = sig[1, 0]
                            stress_component = sig_vect[g_dir]
                        elif self.field.strain_type == StrainType.DISPLACEMENT_SYMMETRIC_GRADIENT:
                            stress_component = material.mat_data.s1.thermodynamic_forces[qp, g_dir]
                        res_qdp_file.write("{},".format(stress_component))
                    qp += 1
                    res_qdp_file.write("\n")
        return

    def solve_newton_1(self, material: Material):
        """

        Args:
            material:
            reset_displacement_at_time_step:

        Returns:

        """
        self.clean_res_dir()
        _dx = self.field.field_dimension
        _fk = self.finite_element.face_basis_k.dimension
        _cl = self.finite_element.cell_basis_l.dimension
        external_forces_coefficient = 1.0
        normalization_lagrange_coefficient = material.lagrange_parameter
        # ----------------------------------------------------------------------------------------------------------
        # SET SYSTEM SIZE
        # ----------------------------------------------------------------------------------------------------------
        _constrained_system_size, _system_size = systm.get_total_system_size(
            self.field, self.finite_element, self.mesh, self.boundary_conditions
        )
        unknown_vector = np.zeros((_constrained_system_size))
        for time_step_index, time_step in enumerate(self.time_steps):
            # --- SET TEMPERATURE
            material.set_temperature()
            correction = np.zeros((_constrained_system_size))
            # --- PRINT DATA
            print("||||||||||||||||||||||||||||||||||")
            print("TIME_STEP : {}".format(time_step))
            # --- WRITE RES FILES
            residual_values = []
            file_suffix = "{}".format(time_step_index).zfill(6)
            self.create_vertex_res_files(file_suffix)
            self.create_quadrature_points_res_files(file_suffix)
            for iteration in range(self.number_of_iterations):
                print("=======================")
                print("ITERATION : {}".format(iteration))
                # --------------------------------------------------------------------------------------------------
                # SET SYSTEM MATRIX AND VECTOR
                # --------------------------------------------------------------------------------------------------
                tangent_matrix = np.zeros((_constrained_system_size, _constrained_system_size))
                residual = np.zeros((_constrained_system_size))
                # --------------------------------------------------------------------------------------------------
                # SET TIME INCREMENT
                # --------------------------------------------------------------------------------------------------
                if time_step_index == 0:
                    _dt = time_step
                else:
                    _dt = time_step - self.time_steps[time_step_index - 1]
                _dt = 0.0
                # --------------------------------------------------------------------------------------------------
                # FOR ELEMENT LOOP
                # --------------------------------------------------------------------------------------------------
                _qp = 0
                for _element_index, element in enumerate(self.elements):
                    # --- GET ELEMENT UNKNOWN CORRECTION
                    # element_unknown_correction = np.copy(element.get_element_unknwon_correction(
                    #     # self.field, self.finite_element, unknown_vector
                    #     self.field, self.finite_element, correction
                    # ))
                    element_unknown_correction = element.get_element_unknwon_correction(
                        self.field, self.finite_element, correction
                    )
                    element.element_unknown_vector += element_unknown_correction
                    print("-------------")
                    print("ELEM : {}".format(_element_index))
                    # print(" * ELEM_CORRECTION : \n{}".format(element_unknown_correction))
                    # print(" * ELEM_CURRENT_DISPLACEMENT : \n{}".format(element.element_unknown_vector))
                    # --- INITIALIZE MATRIX AND VECTORS
                    element_stiffness_matrix = np.zeros((element.element_size, element.element_size))
                    element_internal_forces = np.zeros((element.element_size,))
                    element_external_forces = np.zeros((element.element_size,))
                    # --- RUN OVER EACH QUADRATURE POINT
                    for _qc in range(len(element.cell.quadrature_weights)):
                        _w_q_c = element.cell.quadrature_weights[_qc]
                        _x_q_c = element.cell.quadrature_points[:, _qc]
                        u_at_gauss = np.zeros((_dx,))
                        v = self.finite_element.cell_basis_l.evaluate_function(
                            _x_q_c, element.cell.shape.centroid, element.cell.shape.diameter
                        )
                        for eucli_dir in range(_dx):
                            row_0 = eucli_dir * _cl
                            row_1 = (eucli_dir + 1) * _cl
                            u_at_gauss += v @ element.element_unknown_vector[row_0: row_1]
                        print("U_AT_GAUSS {} : {}".format(_qp, u_at_gauss))
                        # --- COMPUTE STRAINS AND SET THEM IN THE BEHAVIOUR LAW
                        transformation_gradient = element.compute_transformation_gradient(
                            # self.field, _qc, element_unknown_increment
                            self.field, _qc, element.element_unknown_vector
                        )
                        material.mat_data.s1.gradients[_qp] = transformation_gradient
                        # --- INTEGRATE BEHAVIOUR LAW
                        integ_res = mgis_bv.integrate(material.mat_data, material.integration_type, _dt, _qp, (_qp + 1))
                        # --- VOLUMETRIC FORCES
                        v = self.finite_element.cell_basis_l.evaluate_function(
                            _x_q_c, element.cell.shape.centroid, element.cell.shape.diameter
                        )
                        for load in self.loads:
                            vl = _w_q_c * v * load.function(time_step, _x_q_c)
                            _re0 = load.direction * _cl
                            _re1 = (load.direction + 1) * _cl
                            element_external_forces[_re0:_re1] += vl
                        # --- COMPUTE STIFFNESS MATRIX CONTRIBUTION AT QUADRATURE POINT
                        element_stiffness_matrix += _w_q_c * (
                            element.gradients_operators[_qc].T
                            @ material.mat_data.K[_qp]
                            @ element.gradients_operators[_qc]
                        )
                        # --- COMPUTE STIFFNESS MATRIX CONTRIBUTION AT QUADRATURE POINT
                        element_internal_forces += _w_q_c * (
                            element.gradients_operators[_qc].T @ material.mat_data.s1.thermodynamic_forces[_qp]
                        )
                        _qp += 1
                    # --- ADDING STABILIZATION CONTRIBUTION AT THE ELEMENT LEVEL
                    element_stiffness_matrix += material.stabilization_parameter * element.stabilization_operator
                    # --- ADDING STABILIZATION CONTRIBUTION AT THE ELEMENT LEVEL
                    element_internal_forces += (
                        material.stabilization_parameter * element.stabilization_operator @ element.element_unknown_vector
                    )
                    # --- COMPUTE RESIDUAL AFTER VOLUMETRIC CONTRIBUTION
                    element_residual = element_internal_forces - element_external_forces
                    # --- CONDENSATION
                    K_cond, R_cond = element.make_condensation(
                        self.field, self.finite_element, element_stiffness_matrix, element_residual
                    )
                    # --- ASSEMBLY
                    for _i_local, _i_global in enumerate(element.faces_indices):
                        _rg0 = (_dx * _fk) * _i_global
                        _rg1 = (_dx * _fk) * (_i_global + 1)
                        _re0 = _i_local * (_dx * _fk)
                        _re1 = (_i_local + 1) * (_dx * _fk)
                        residual[_rg0:_rg1] += R_cond[_re0:_re1]
                        for _j_local, _j_global in enumerate(element.faces_indices):
                            _cg0 = _j_global * (_dx * _fk)
                            _cg1 = (_j_global + 1) * (_dx * _fk)
                            _ce0 = _j_local * (_dx * _fk)
                            _ce1 = (_j_local + 1) * (_dx * _fk)
                            tangent_matrix[_rg0:_rg1, _cg0:_cg1] += K_cond[_re0:_re1, _ce0:_ce1]
                    # --- SET EXTERNAL FORCES COEFFICIENT
                    if np.max(np.abs(element_external_forces)) > external_forces_coefficient:
                        external_forces_coefficient = np.max(np.abs(element_external_forces))
                # --------------------------------------------------------------------------------------------------
                # FOR FACE LOOP
                # --------------------------------------------------------------------------------------------------
                iter_face_constraint = 0
                for boundary_condition in self.boundary_conditions:
                    # --- DISPLACEMENT CONDITIONS
                    if boundary_condition.boundary_type == BoundaryType.DISPLACEMENT:
                        for element in self.elements:
                            # element_unknown_increment = element.get_element_unknwon_increment(
                            #     self.field, self.finite_element, unknown_increment
                            # )
                            for f_local, f_global in enumerate(element.faces_indices):
                                if (
                                    f_global
                                    in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]
                                ):
                                    _rl0 = _system_size + iter_face_constraint * _fk
                                    _rl1 = _system_size + (iter_face_constraint + 1) * _fk
                                    _re0 = _cl * _dx + (f_local * _dx * _fk) + boundary_condition.direction * _fk
                                    _re1 = _cl * _dx + (f_local * _dx * _fk) + (boundary_condition.direction + 1) * _fk
                                    _rg0 = f_global * _fk * _dx + _fk * boundary_condition.direction
                                    _rg1 = f_global * _fk * _dx + _fk * (boundary_condition.direction + 1)
                                    # -------
                                    # face_displacement = element_unknown_increment[_r0:_r1]
                                    # face_displacement = np.copy(element_unknown_increment[_c0:_c1])
                                    # face_displacement = element_unknown_increment[_c0:_c1]
                                    face_lagrange = unknown_vector[_rl0:_rl1]
                                    face_displacement = unknown_vector[_rg0:_rg1]
                                    _m_psi_psi_face = np.zeros((_fk, _fk), dtype=real)
                                    _v_face_imposed_displacement = np.zeros((_fk,), dtype=real)
                                    for _qf in range(len(element.faces[f_local].quadrature_weights)):
                                        _x_q_f = element.faces[f_local].quadrature_points[:, _qf]
                                        _w_q_f = element.faces[f_local].quadrature_weights[_qf]
                                        v = self.finite_element.face_basis_k.evaluate_function(
                                            _x_q_f,
                                            element.faces[f_local].shape.centroid,
                                            element.faces[f_local].shape.diameter,
                                        )
                                        _v_face_imposed_displacement += (
                                            _w_q_f * v * boundary_condition.function(time_step, _x_q_f)
                                        )
                                        _m_psi_psi_face += blocks.get_face_mass_matrix_in_face(
                                            element.faces[f_local],
                                            self.finite_element.face_basis_k,
                                            self.finite_element.face_basis_k,
                                            _x_q_f,
                                            _w_q_f,
                                        )
                                    _m_psi_psi_face_inv = np.linalg.inv(_m_psi_psi_face)
                                    imposed_face_displacement = _m_psi_psi_face_inv @ _v_face_imposed_displacement
                                    face_displacement_difference = face_displacement - imposed_face_displacement
                                    # --- LAGRANGE INTERNAL FORCES PART
                                    residual[_rg0:_rg1] += normalization_lagrange_coefficient * face_lagrange
                                    # --- LAGRANGE MULTIPLIERS PART
                                    residual[_rl0:_rl1] += (
                                        normalization_lagrange_coefficient * face_displacement_difference
                                    )
                                    # --- LAGRANGE MATRIX PART
                                    tangent_matrix[_rl0:_rl1, _rg0:_rg1] += normalization_lagrange_coefficient * np.eye(
                                        _fk
                                    )
                                    tangent_matrix[_rg0:_rg1, _rl0:_rl1] += normalization_lagrange_coefficient * np.eye(
                                        _fk
                                    )
                                    # --- SET EXTERNAL FORCES COEFFICIENT
                                    lagrange_external_forces = (
                                        normalization_lagrange_coefficient * imposed_face_displacement
                                    )
                                    if np.max(np.abs(lagrange_external_forces)) > external_forces_coefficient:
                                        external_forces_coefficient = np.max(np.abs(lagrange_external_forces))
                                    # if np.max(np.abs(residual[_rg0:_rg1])) > external_forces_coefficient:
                                    #     external_forces_coefficient = np.max(np.abs(residual[_rg0:_rg1]))
                                    # --- ITER ON LAGRANGE CONSTRAINED FACE
                                    iter_face_constraint += 1
                    elif boundary_condition.boundary_type == BoundaryType.PRESSURE:
                        for element in self.elements:
                            # element_unknown_increment = element.get_element_unknwon_increment(
                            #     self.field, self.finite_element, unknown_increment
                            # )
                            for f_local, f_global in enumerate(element.faces_indices):
                                if (
                                    f_global
                                    in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]
                                ):
                                    for qf in range(len(element.faces[f_local].quadrature_weights)):
                                        _x_q_f = element.faces[f_local].quadrature_points[:, qf]
                                        _w_q_f = element.faces[f_local].quadrature_weights[qf]
                                        v = self.finite_element.face_basis_k.evaluate_function(
                                            _x_q_f,
                                            element.faces[f_local].shape.centroid,
                                            element.faces[f_local].shape.diameter,
                                        )
                                        vf = _w_q_f * v * boundary_condition.function(time_step, _x_q_f)
                                        # _c0 = _dx * _cl + f_local * _dx * _fk + boundary_condition.direction * _fk
                                        # _c1 = _dx * _cl + f_local * _dx * _fk + (boundary_condition.direction + 1) * _fk
                                        _rg0 = f_global * _fk * _dx + _fk * boundary_condition.direction
                                        _rg1 = f_global * _fk * _dx + _fk * (boundary_condition.direction + 1)
                                        # residual[_c0:_c1] += vf
                                        residual[_rg0:_rg1] -= vf
                                        # --- SET EXTERNAL FORCES COEFFICIENT
                                        if np.max(np.abs(vf)) > external_forces_coefficient:
                                            external_forces_coefficient = np.max(np.abs(vf))
                # --------------------------------------------------------------------------------------------------
                # RESIDUAL EVALUATION
                # --------------------------------------------------------------------------------------------------
                tol_vect = np.ones((_constrained_system_size,), dtype=real) * self.tolerance
                if external_forces_coefficient == 0.0:
                    external_forces_coefficient = 1.0
                res_eval = (1.0 / external_forces_coefficient) * residual
                print("ITER : {} | RES_MAX : {}".format(str(iteration).zfill(4), max(np.abs(res_eval))))
                # print("------------")
                residual_values.append(max(np.abs(res_eval)))
                if (np.abs(res_eval) < tol_vect).all():
                    # ----------------------------------------------------------------------------------------------
                    # UPDATE INTERNAL VARIABLES
                    # ----------------------------------------------------------------------------------------------
                    mgis_bv.update(material.mat_data)
                    # print("ITERATIONS : {}".format(iteration + 1))
                    # print("|||||||||||||||||||||")
                    self.write_vertex_res_files(file_suffix, unknown_vector)
                    self.write_quadrature_points_res_files(file_suffix, unknown_vector, material)
                    # plt.plot(range(len(residual_values)), residual_values)
                    # plt.show()
                    break
                else:
                    # ----------------------------------------------------------------------------------------------
                    # SOLVE SYSTEM
                    # ----------------------------------------------------------------------------------------------
                    if iteration == self.number_of_iterations - 1:
                        return
                    else:
                        sparse_global_matrix = csr_matrix(-tangent_matrix)
                        correction = spsolve(sparse_global_matrix, residual)
                        unknown_vector += correction

        return

    def solve_newton_exact(self, material: Material, reset_displacement_at_time_step: bool):
        """

        Args:
            material:
            reset_displacement_at_time_step:

        Returns:

        """
        self.clean_res_dir()
        _dx = self.field.field_dimension
        _fk = self.finite_element.face_basis_k.dimension
        _cl = self.finite_element.cell_basis_l.dimension
        external_forces_coefficient = 1.0
        normalization_lagrange_coefficient = material.lagrange_parameter
        # ----------------------------------------------------------------------------------------------------------
        # SET SYSTEM SIZE
        # ----------------------------------------------------------------------------------------------------------
        _constrained_system_size, _system_size = systm.get_total_system_size(
            self.field, self.finite_element, self.mesh, self.boundary_conditions
        )
        unknown_increment = np.zeros((_constrained_system_size))
        for time_step_index, time_step in enumerate(self.time_steps):
            # --- SET TEMPERATURE
            if reset_displacement_at_time_step:
                unknown_increment = np.zeros((_constrained_system_size))
            material.set_temperature()
            # --- PRINT DATA
            print("-------------")
            print("TIME_STEP : {}".format(time_step))
            # --- WRITE RES FILES
            residual_values = []
            file_suffix = "{}".format(time_step_index).zfill(6)
            self.create_vertex_res_files(file_suffix)
            self.create_quadrature_points_res_files(file_suffix)
            for iteration in range(self.number_of_iterations):
                # --------------------------------------------------------------------------------------------------
                # SET SYSTEM MATRIX AND VECTOR
                # --------------------------------------------------------------------------------------------------
                tangent_matrix = np.zeros((_constrained_system_size, _constrained_system_size))
                residual = np.zeros((_constrained_system_size))
                # --------------------------------------------------------------------------------------------------
                # SET TIME INCREMENT
                # --------------------------------------------------------------------------------------------------
                if time_step_index == 0:
                    _dt = time_step
                else:
                    _dt = time_step - self.time_steps[time_step_index - 1]
                _dt = 0.0
                # --------------------------------------------------------------------------------------------------
                # FOR ELEMENT LOOP
                # --------------------------------------------------------------------------------------------------
                _qp = 0
                for element in self.elements:
                    # --- GET ELEMENT UNKNOWN INCREMENT
                    element_unknown_increment = element.get_element_unknwon_correction(
                        self.field, self.finite_element, unknown_increment
                    )
                    # --- INITIALIZE MATRIX AND VECTORS
                    element_stiffness_matrix = np.zeros((element.element_size, element.element_size))
                    element_internal_forces = np.zeros((element.element_size,))
                    element_external_forces = np.zeros((element.element_size,))
                    _cs = _cl * _dx
                    # _nf = (len(element.faces_indices))
                    # _fs = _nf * _fk * _dx
                    # cell_unknown_increment = np.zeros((_cs,))
                    # face_unknown_increment = np.zeros((_fs,))
                    # face_unknown_increment = element.get_face_unknown_increment(self.field, self.finite_element, unknown_increment)
                    local_element_unknown_increment = np.copy(element_unknown_increment)
                    cell_iterations = 10
                    for cell_iteration in range(cell_iterations):
                        # --- RUN OVER EACH QUADRATURE POINT
                        _qp_count = 0
                        for _qc in range(len(element.cell.quadrature_weights)):
                            _w_q_c = element.cell.quadrature_weights[_qc]
                            _x_q_c = element.cell.quadrature_points[:, _qc]
                            # --- COMPUTE STRAINS AND SET THEM IN THE BEHAVIOUR LAW
                            # local_element_unknown_increment
                            # transformation_gradient = element.compute_transformation_gradient(
                            #     self.field, _qc, element_unknown_increment
                            # )
                            transformation_gradient = element.compute_transformation_gradient(
                                self.field, _qc, local_element_unknown_increment
                            )
                            material.mat_data.s1.gradients[_qp] = transformation_gradient
                            # --- INTEGRATE BEHAVIOUR LAW
                            integ_res = mgis_bv.integrate(material.mat_data, material.integration_type, _dt, _qp, (_qp + 1))
                            # --- VOLUMETRIC FORCES
                            v = self.finite_element.cell_basis_l.evaluate_function(
                                _x_q_c, element.cell.shape.centroid, element.cell.shape.diameter
                            )
                            for load in self.loads:
                                vl = _w_q_c * v * load.function(time_step, _x_q_c)
                                _re0 = load.direction * _cl
                                _re1 = (load.direction + 1) * _cl
                                element_external_forces[_re0:_re1] += vl
                            # --- COMPUTE STIFFNESS MATRIX CONTRIBUTION AT QUADRATURE POINT
                            element_stiffness_matrix += _w_q_c * (
                                element.gradients_operators[_qc].T
                                @ material.mat_data.K[_qp]
                                @ element.gradients_operators[_qc]
                            )
                            # --- COMPUTE STIFFNESS MATRIX CONTRIBUTION AT QUADRATURE POINT
                            element_internal_forces += _w_q_c * (
                                element.gradients_operators[_qc].T @ material.mat_data.s1.thermodynamic_forces[_qp]
                            )
                            _qp += 1
                            _qp_count += 1
                        # if cell_iteration == cell_iterations - 1:
                        _qp -= _qp_count
                        # --- ADDING STABILIZATION CONTRIBUTION AT THE ELEMENT LEVEL
                        element_stiffness_matrix += material.stabilization_parameter * element.stabilization_operator
                        # --- ADDING STABILIZATION CONTRIBUTION AT THE ELEMENT LEVEL
                        element_internal_forces += (
                            material.stabilization_parameter * element.stabilization_operator @ element_unknown_increment
                        )
                        # --- COMPUTE RESIDUAL AFTER VOLUMETRIC CONTRIBUTION
                        element_residual = element_internal_forces - element_external_forces
                        # --- HERE
                        K_cc = element_stiffness_matrix[:_cs, :_cs]
                        K_cf = element_stiffness_matrix[:_cs, _cs:]
                        K_fc = element_stiffness_matrix[_cs:, :_cs]
                        K_ff = element_stiffness_matrix[_cs:, _cs:]
                        K_cc_inv = np.linalg.inv(K_cc)
                        R_c = element_residual[:_cs]
                        # cell_unknown_increment = np.linalg.solve(-K_cc, R_c)
                        cell_unknown_correction = -K_cc_inv @ R_c
                        local_element_unknown_increment[:_cs] += cell_unknown_correction
                        print("LOCAL_CELL_RES : {}".format(np.max(np.abs(R_c))))
                    # --- END CELL NEWTON
                    # local_element_unknown_increment[_cs:] = (K_ff - K_fc @ K_cc_inv @ K_cf) @ element_residual[_cs:]
                    face_stiffness_matrix = (K_ff - K_fc @ K_cc_inv @ K_cf)
                    face_residual = element_residual[_cs:]
                    print("LOCAL_FACE_RES : {}".format(np.max(np.abs(face_residual))))
                    # local_face_unknwon_increment = (K_ff - K_fc @ K_cc_inv @ K_cf) @ element_residual[_cs:]
                    for _i_local, _i_global in enumerate(element.faces_indices):
                        _rg0 = (_dx * _fk) * _i_global
                        _rg1 = (_dx * _fk) * (_i_global + 1)
                        _re0 = _i_local * (_dx * _fk)
                        _re1 = (_i_local + 1) * (_dx * _fk)
                        residual[_rg0:_rg1] += face_residual[_re0:_re1]
                        for _j_local, _j_global in enumerate(element.faces_indices):
                            _cg0 = _j_global * (_dx * _fk)
                            _cg1 = (_j_global + 1) * (_dx * _fk)
                            _ce0 = _j_local * (_dx * _fk)
                            _ce1 = (_j_local + 1) * (_dx * _fk)
                            tangent_matrix[_rg0:_rg1, _cg0:_cg1] += face_stiffness_matrix[_re0:_re1, _ce0:_ce1]
                    # --- WRITE FACE INCREMENT ON GLOBAL UNKNOWN INCREMENT


                    # # --- CONDENSATION
                    # K_cond, R_cond = element.make_condensation(
                    #     self.field, self.finite_element, element_stiffness_matrix, element_residual
                    # )
                    # # --- ASSEMBLY
                    # for _i_local, _i_global in enumerate(element.faces_indices):
                    #     _rg0 = (_dx * _fk) * _i_global
                    #     _rg1 = (_dx * _fk) * (_i_global + 1)
                    #     _re0 = _i_local * (_dx * _fk)
                    #     _re1 = (_i_local + 1) * (_dx * _fk)
                    #     residual[_rg0:_rg1] += R_cond[_re0:_re1]
                    #     for _j_local, _j_global in enumerate(element.faces_indices):
                    #         _cg0 = _j_global * (_dx * _fk)
                    #         _cg1 = (_j_global + 1) * (_dx * _fk)
                    #         _ce0 = _j_local * (_dx * _fk)
                    #         _ce1 = (_j_local + 1) * (_dx * _fk)
                    #         tangent_matrix[_rg0:_rg1, _cg0:_cg1] += K_cond[_re0:_re1, _ce0:_ce1]
                    # --- SET EXTERNAL FORCES COEFFICIENT
                    if np.max(np.abs(element_external_forces)) > external_forces_coefficient:
                        external_forces_coefficient = np.max(np.abs(element_external_forces))
                # --------------------------------------------------------------------------------------------------
                # FOR FACE LOOP
                # --------------------------------------------------------------------------------------------------
                iter_face_constraint = 0
                for boundary_condition in self.boundary_conditions:
                    # --- DISPLACEMENT CONDITIONS
                    if boundary_condition.boundary_type == BoundaryType.DISPLACEMENT:
                        for element in self.elements:
                            # element_unknown_increment = element.get_element_unknwon_increment(
                            #     self.field, self.finite_element, unknown_increment
                            # )
                            for f_local, f_global in enumerate(element.faces_indices):
                                if (
                                    f_global
                                    in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]
                                ):
                                    _rl0 = _system_size + iter_face_constraint * _fk
                                    _rl1 = _system_size + (iter_face_constraint + 1) * _fk
                                    _re0 = _cl * _dx + (f_local * _dx * _fk) + boundary_condition.direction * _fk
                                    _re1 = _cl * _dx + (f_local * _dx * _fk) + (boundary_condition.direction + 1) * _fk
                                    _rg0 = f_global * _fk * _dx + _fk * boundary_condition.direction
                                    _rg1 = f_global * _fk * _dx + _fk * (boundary_condition.direction + 1)
                                    # -------
                                    # face_displacement = element_unknown_increment[_r0:_r1]
                                    # face_displacement = np.copy(element_unknown_increment[_c0:_c1])
                                    # face_displacement = element_unknown_increment[_c0:_c1]
                                    face_lagrange = unknown_increment[_rl0:_rl1]
                                    face_displacement = unknown_increment[_rg0:_rg1]
                                    _m_psi_psi_face = np.zeros((_fk, _fk), dtype=real)
                                    _v_face_imposed_displacement = np.zeros((_fk,), dtype=real)
                                    for _qf in range(len(element.faces[f_local].quadrature_weights)):
                                        _x_q_f = element.faces[f_local].quadrature_points[:, _qf]
                                        _w_q_f = element.faces[f_local].quadrature_weights[_qf]
                                        v = self.finite_element.face_basis_k.evaluate_function(
                                            _x_q_f,
                                            element.faces[f_local].shape.centroid,
                                            element.faces[f_local].shape.diameter,
                                        )
                                        _v_face_imposed_displacement += (
                                            _w_q_f * v * boundary_condition.function(time_step, _x_q_f)
                                        )
                                        _m_psi_psi_face += blocks.get_face_mass_matrix_in_face(
                                            element.faces[f_local],
                                            self.finite_element.face_basis_k,
                                            self.finite_element.face_basis_k,
                                            _x_q_f,
                                            _w_q_f,
                                        )
                                    _m_psi_psi_face_inv = np.linalg.inv(_m_psi_psi_face)
                                    imposed_face_displacement = _m_psi_psi_face_inv @ _v_face_imposed_displacement
                                    face_displacement_difference = face_displacement - imposed_face_displacement
                                    # --- LAGRANGE INTERNAL FORCES PART
                                    residual[_rg0:_rg1] += normalization_lagrange_coefficient * face_lagrange
                                    # --- LAGRANGE MULTIPLIERS PART
                                    residual[_rl0:_rl1] += (
                                        normalization_lagrange_coefficient * face_displacement_difference
                                    )
                                    # --- LAGRANGE MATRIX PART
                                    tangent_matrix[_rl0:_rl1, _rg0:_rg1] += normalization_lagrange_coefficient * np.eye(
                                        _fk
                                    )
                                    tangent_matrix[_rg0:_rg1, _rl0:_rl1] += normalization_lagrange_coefficient * np.eye(
                                        _fk
                                    )
                                    # --- SET EXTERNAL FORCES COEFFICIENT
                                    lagrange_external_forces = (
                                        normalization_lagrange_coefficient * imposed_face_displacement
                                    )
                                    if np.max(np.abs(lagrange_external_forces)) > external_forces_coefficient:
                                        external_forces_coefficient = np.max(np.abs(lagrange_external_forces))
                                    # if np.max(np.abs(residual[_rg0:_rg1])) > external_forces_coefficient:
                                    #     external_forces_coefficient = np.max(np.abs(residual[_rg0:_rg1]))
                                    # --- ITER ON LAGRANGE CONSTRAINED FACE
                                    iter_face_constraint += 1
                    elif boundary_condition.boundary_type == BoundaryType.PRESSURE:
                        for element in self.elements:
                            # element_unknown_increment = element.get_element_unknwon_increment(
                            #     self.field, self.finite_element, unknown_increment
                            # )
                            for f_local, f_global in enumerate(element.faces_indices):
                                if (
                                    f_global
                                    in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]
                                ):
                                    for qf in range(len(element.faces[f_local].quadrature_weights)):
                                        _x_q_f = element.faces[f_local].quadrature_points[:, qf]
                                        _w_q_f = element.faces[f_local].quadrature_weights[qf]
                                        v = self.finite_element.face_basis_k.evaluate_function(
                                            _x_q_f,
                                            element.faces[f_local].shape.centroid,
                                            element.faces[f_local].shape.diameter,
                                        )
                                        vf = _w_q_f * v * boundary_condition.function(time_step, _x_q_f)
                                        # _c0 = _dx * _cl + f_local * _dx * _fk + boundary_condition.direction * _fk
                                        # _c1 = _dx * _cl + f_local * _dx * _fk + (boundary_condition.direction + 1) * _fk
                                        _rg0 = f_global * _fk * _dx + _fk * boundary_condition.direction
                                        _rg1 = f_global * _fk * _dx + _fk * (boundary_condition.direction + 1)
                                        # residual[_c0:_c1] += vf
                                        residual[_rg0:_rg1] -= vf
                                        # --- SET EXTERNAL FORCES COEFFICIENT
                                        if np.max(np.abs(vf)) > external_forces_coefficient:
                                            external_forces_coefficient = np.max(np.abs(vf))
                # --------------------------------------------------------------------------------------------------
                # RESIDUAL EVALUATION
                # --------------------------------------------------------------------------------------------------
                tol_vect = np.ones((_constrained_system_size,), dtype=real) * self.tolerance
                if external_forces_coefficient == 0.0:
                    external_forces_coefficient = 1.0
                res_eval = (1.0 / external_forces_coefficient) * residual
                print("ITER : {} | RES_MAX : {}".format(str(iteration).zfill(4), max(np.abs(res_eval))))
                residual_values.append(max(np.abs(res_eval)))
                if (np.abs(res_eval) < tol_vect).all():
                    # ----------------------------------------------------------------------------------------------
                    # UPDATE INTERNAL VARIABLES
                    # ----------------------------------------------------------------------------------------------
                    mgis_bv.update(material.mat_data)
                    print("ITERATIONS : {}".format(iteration + 1))
                    self.write_vertex_res_files(file_suffix, unknown_increment)
                    self.write_quadrature_points_res_files(file_suffix, unknown_increment, material)
                    # plt.plot(range(len(residual_values)), residual_values)
                    # plt.show()
                    break
                else:
                    # ----------------------------------------------------------------------------------------------
                    # SOLVE SYSTEM
                    # ----------------------------------------------------------------------------------------------
                    sparse_global_matrix = csr_matrix(-tangent_matrix)
                    correction = spsolve(sparse_global_matrix, residual)
                    unknown_increment += correction
        return

    def solve_newton_check_1(self, material: Material, reset_displacement_at_time_step: bool):
        """

        Args:
            material:
            reset_displacement_at_time_step:

        Returns:

        """
        self.clean_res_dir()
        _dx = self.field.field_dimension
        _fk = self.finite_element.face_basis_k.dimension
        _cl = self.finite_element.cell_basis_l.dimension
        external_forces_coefficient = 1.0
        normalization_lagrange_coefficient = material.lagrange_parameter
        # ----------------------------------------------------------------------------------------------------------
        # SET SYSTEM SIZE
        # ----------------------------------------------------------------------------------------------------------
        _constrained_system_size, _system_size = systm.get_total_system_size(
            self.field, self.finite_element, self.mesh, self.boundary_conditions
        )
        unknown_increment = np.zeros((_constrained_system_size))
        """correction = np.zeros((_constrained_system_size))
        """
        for time_step_index, time_step in enumerate(self.time_steps):
            # --- SET TEMPERATURE
            if reset_displacement_at_time_step:
                unknown_increment = np.zeros((_constrained_system_size))
            material.set_temperature()
            # --- PRINT DATA
            print("-------------")
            print("TIME_STEP : {}".format(time_step))
            # --- WRITE RES FILES
            residual_values = []
            file_suffix = "{}".format(time_step_index).zfill(6)
            self.create_vertex_res_files(file_suffix)
            self.create_quadrature_points_res_files(file_suffix)
            for iteration in range(self.number_of_iterations):
                # --------------------------------------------------------------------------------------------------
                # SET SYSTEM MATRIX AND VECTOR
                # --------------------------------------------------------------------------------------------------
                tangent_matrix = np.zeros((_constrained_system_size, _constrained_system_size))
                residual = np.zeros((_constrained_system_size))
                # --------------------------------------------------------------------------------------------------
                # SET TIME INCREMENT
                # --------------------------------------------------------------------------------------------------
                if time_step_index == 0:
                    _dt = time_step
                else:
                    _dt = time_step - self.time_steps[time_step_index - 1]
                _dt = 0.0
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                noise_tolerance = 1.0e-6
                noise = noise_tolerance * np.max(np.abs(np.copy(unknown_increment[:_system_size])))
                # if noise == 0.:
                #     noise = 1.e9
                numerical_forces_bs_0 = [np.zeros((element.element_size, element.element_size), dtype=real) for element in self.elements]
                numerical_forces_as_0 = [np.zeros((element.element_size, element.element_size), dtype=real) for element in self.elements]
                numerical_unknowns_0 = [np.zeros((element.element_size, element.element_size), dtype=real) for element in self.elements]
                for epsi_dir in range(self.elements[0].element_size):
                    _qp = 0
                    for _element_index, element in enumerate(self.elements):
                        # --- GET ELEMENT UNKNOWN INCREMENT
                        element_unknown_increment = element.get_element_unknwon_correction(
                            self.field, self.finite_element, unknown_increment
                        )
                        numerical_unknowns_0[_element_index][:, epsi_dir] += np.copy(element_unknown_increment)
                        numerical_unknowns_0[_element_index][epsi_dir, epsi_dir] += noise
                        for _qc in range(len(element.cell.quadrature_weights)):
                            _w_q_c = element.cell.quadrature_weights[_qc]
                            _x_q_c = element.cell.quadrature_points[:, _qc]
                            transformation_gradient_0 = element.compute_transformation_gradient(
                                self.field, _qc, numerical_unknowns_0[_element_index][:, epsi_dir]
                            )
                            material.mat_data.s1.gradients[_qp] = transformation_gradient_0
                            integ_res = mgis_bv.integrate(
                                material.mat_data, material.integration_type, _dt, _qp, (_qp + 1)
                            )
                            numerical_forces_bs_0[_element_index][:, epsi_dir] += _w_q_c * (
                                element.gradients_operators[_qc].T @ material.mat_data.s1.thermodynamic_forces[_qp]
                            )
                            _qp += 1
                                # mgis_bv.revert(material.mat_data)
                        # numerical_matrices_0_bs.append(np.copy(element_internal_epsi_0))
                    mgis_bv.revert(material.mat_data)
                for _element_index, element in enumerate(self.elements):
                    numerical_forces_as_0[_element_index] += numerical_forces_bs_0[_element_index]
                    for epsi_dir in range(self.elements[0].element_size):
                        numerical_forces_as_0[_element_index][:, epsi_dir] += (
                                material.stabilization_parameter
                                * element.stabilization_operator
                                @ numerical_unknowns_0[_element_index][:, epsi_dir]
                        )
                    # for _element_index, element in enumerate(self.elements):
                    #     numerical_forces_as_0[_element_index] = numerical_forces_bs_0[_element_index]
                    #     numerical_forces_as_0[_element_index][:, epsi_dir] += (
                    #             material.stabilization_parameter
                    #             * element.stabilization_operator
                    #             @ numerical_unknowns_0[_element_index][:, epsi_dir]
                    #     )
                numerical_forces_bs_1 = [np.zeros((element.element_size, element.element_size), dtype=real) for element
                                         in self.elements]
                numerical_forces_as_1 = [np.zeros((element.element_size, element.element_size), dtype=real) for element
                                         in self.elements]
                numerical_unknowns_1 = [np.zeros((element.element_size, element.element_size), dtype=real) for element
                                        in self.elements]
                for epsi_dir in range(self.elements[0].element_size):
                    _qp = 0
                    for _element_index, element in enumerate(self.elements):
                        # --- GET ELEMENT UNKNOWN INCREMENT
                        element_unknown_increment = element.get_element_unknwon_correction(
                            self.field, self.finite_element, unknown_increment
                        )
                        numerical_unknowns_1[_element_index][:, epsi_dir] += np.copy(element_unknown_increment)
                        numerical_unknowns_1[_element_index][epsi_dir, epsi_dir] -= noise
                        for _qc in range(len(element.cell.quadrature_weights)):
                            _w_q_c = element.cell.quadrature_weights[_qc]
                            _x_q_c = element.cell.quadrature_points[:, _qc]
                            transformation_gradient_0 = element.compute_transformation_gradient(
                                self.field, _qc, numerical_unknowns_1[_element_index][:, epsi_dir]
                            )
                            material.mat_data.s1.gradients[_qp] = transformation_gradient_0
                            integ_res = mgis_bv.integrate(
                                material.mat_data, material.integration_type, _dt, _qp, (_qp + 1)
                            )
                            numerical_forces_bs_1[_element_index][:, epsi_dir] += _w_q_c * (
                                    element.gradients_operators[_qc].T @ material.mat_data.s1.thermodynamic_forces[_qp]
                            )
                            _qp += 1
                            # mgis_bv.revert(material.mat_data)
                        # numerical_matrices_0_bs.append(np.copy(element_internal_epsi_0))
                    mgis_bv.revert(material.mat_data)
                    # for _element_index, element in enumerate(self.elements):
                    #     numerical_forces_as_1[_element_index] += numerical_forces_bs_1[_element_index]
                    #     numerical_forces_as_1[_element_index][:, epsi_dir] += (
                    #             material.stabilization_parameter
                    #             * element.stabilization_operator
                    #             @ numerical_unknowns_1[_element_index][:, epsi_dir]
                    #     )
                for _element_index, element in enumerate(self.elements):
                    numerical_forces_as_1[_element_index] = +numerical_forces_bs_1[_element_index]
                    for epsi_dir in range(self.elements[0].element_size):
                        numerical_forces_as_1[_element_index][:, epsi_dir] += (
                                material.stabilization_parameter
                                * element.stabilization_operator
                                @ numerical_unknowns_1[_element_index][:, epsi_dir]
                        )
                # numerical_matrices_1_bs = []
                # numerical_matrices_1_as = []
                # _qp = 0
                # for _element_index, element in enumerate(self.elements):
                #     # --- GET ELEMENT UNKNOWN INCREMENT
                #     element_unknown_increment = element.get_element_unknwon_increment(
                #         self.field, self.finite_element, unknown_increment
                #     )
                #     element_unknown_epsi_1 = np.zeros((element.element_size, element.element_size), dtype=real)
                #     element_internal_epsi_1 = np.zeros((element.element_size, element.element_size), dtype=real)
                #     if noise == 0.:
                #         noise = 1.
                #     for epsi_dir in range(element.element_size):
                #         element_unknown_epsi_1[:, epsi_dir] += np.copy(element_unknown_increment)
                #         element_unknown_epsi_1[epsi_dir, epsi_dir] -= noise
                #     for _qc in range(len(element.cell.quadrature_weights)):
                #         _w_q_c = element.cell.quadrature_weights[_qc]
                #         _x_q_c = element.cell.quadrature_points[:, _qc]
                #         for epsi_dir in range(element.element_size):
                #             transformation_gradient_0 = element.compute_transformation_gradient(
                #                 self.field, _qc, element_unknown_epsi_1[:, epsi_dir]
                #             )
                #             material.mat_data.s1.gradients[_qp] = transformation_gradient_0
                #             integ_res = mgis_bv.integrate(
                #                 material.mat_data, material.integration_type, _dt, _qp, (_qp + 1)
                #             )
                #             element_internal_epsi_1[:, epsi_dir] += _w_q_c * (
                #                     element.gradients_operators[_qc].T @ material.mat_data.s1.thermodynamic_forces[_qp]
                #             )
                #         _qp += 1
                #             # mgis_bv.revert(material.mat_data)
                #     numerical_matrices_1_bs.append(np.copy(element_internal_epsi_1))
                #     # check_threshold = 10.
                #     # numerical_element_stiffness_matrix = (element_internal_epsi_0 - element_internal_epsi_1) / (
                #     #         2.0 * noise
                #     # )
                #     # # check_num = np.max(np.abs(numerical_element_stiffness_matrix)) - np.max(np.abs(element_stiffness_matrix))
                #     # check_num = np.max(np.abs(numerical_element_stiffness_matrix - element_stiffness_matrix))
                #     # if check_num > check_threshold:
                #     #     print("AT ELEM : {} | CHECK_BEFORE_STAB : {}".format(str(_element_index).zfill(6),
                #     #                                                          check_num))
                #     for epsi_dir in range(element.element_size):
                #         element_internal_epsi_1[:, epsi_dir] += (
                #                 material.stabilization_parameter
                #                 * element.stabilization_operator
                #                 @ element_unknown_epsi_1[:, epsi_dir]
                #         )
                #     numerical_matrices_1_as.append(np.copy(element_internal_epsi_1))
                # mgis_bv.revert(material.mat_data)
                    # numerical_element_stiffness_matrix = (element_internal_epsi_0 - element_internal_epsi_1) / (
                    #     2.0 * noise
                    # )
                    # # check_num = np.max(np.abs(numerical_element_stiffness_matrix)) - np.max(np.abs(element_stiffness_matrix))
                    # check_num = np.max(np.abs(numerical_element_stiffness_matrix - element_stiffness_matrix))
                    # if check_num > check_threshold:
                    #     print("AT ELEM : {} | CHECK_AFTER_STAB : {}".format(str(_element_index).zfill(6), check_num))
                # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # --------------------------------------------------------------------------------------------------
                # FOR ELEMENT LOOP
                # --------------------------------------------------------------------------------------------------
                _qp = 0
                for _element_index, element in enumerate(self.elements):
                    # --- GET ELEMENT UNKNOWN INCREMENT
                    element_unknown_increment = element.get_element_unknwon_correction(
                        self.field, self.finite_element, unknown_increment
                    )
                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    # element_unknown_epsi_0 = np.zeros((element.element_size, element.element_size), dtype=real)
                    # element_unknown_epsi_1 = np.zeros((element.element_size, element.element_size), dtype=real)
                    # element_internal_epsi_0 = np.zeros((element.element_size, element.element_size), dtype=real)
                    # element_internal_epsi_1 = np.zeros((element.element_size, element.element_size), dtype=real)
                    # noise = 1.0e-4 * np.max(np.abs(np.copy(unknown_increment[:_system_size])))
                    # if noise == 0.:
                    #     noise = 1.
                    # for epsi_dir in range(element.element_size):
                    #     element_unknown_epsi_0[:, epsi_dir] += np.copy(element_unknown_increment)
                    #     element_unknown_epsi_1[:, epsi_dir] += np.copy(element_unknown_increment)
                    #     element_unknown_epsi_0[epsi_dir, epsi_dir] += noise
                    #     element_unknown_epsi_1[epsi_dir, epsi_dir] -= noise
                    # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # --- INITIALIZE MATRIX AND VECTORS
                    element_stiffness_matrix = np.zeros((element.element_size, element.element_size))
                    element_internal_forces = np.zeros((element.element_size,))
                    element_external_forces = np.zeros((element.element_size,))
                    # --- RUN OVER EACH QUADRATURE POINT
                    for _qc in range(len(element.cell.quadrature_weights)):
                        _w_q_c = element.cell.quadrature_weights[_qc]
                        _x_q_c = element.cell.quadrature_points[:, _qc]
                        # --- COMPUTE STRAINS AND SET THEM IN THE BEHAVIOUR LAW
                        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                        # for epsi_dir in range(element.element_size):
                        #     # --- FIRST PERTURBATION
                        #     transformation_gradient_0 = element.compute_transformation_gradient(
                        #         self.field, _qc, element_unknown_epsi_0[:, epsi_dir]
                        #     )
                        #     material.mat_data.s1.gradients[_qp] = transformation_gradient_0
                        #     integ_res = mgis_bv.integrate(
                        #         material.mat_data, material.integration_type, _dt, _qp, (_qp + 1)
                        #     )
                        #     element_internal_epsi_0[:, epsi_dir] += _w_q_c * (
                        #         element.gradients_operators[_qc].T @ material.mat_data.s1.thermodynamic_forces[_qp]
                        #     )
                        #     mgis_bv.revert(material.mat_data)
                        #     # --- SECOND PERTURBATION
                        #     transformation_gradient_1 = element.compute_transformation_gradient(
                        #         self.field, _qc, element_unknown_epsi_1[:, epsi_dir]
                        #     )
                        #     material.mat_data.s1.gradients[_qp] = transformation_gradient_1
                        #     integ_res = mgis_bv.integrate(
                        #         material.mat_data, material.integration_type, _dt, _qp, (_qp + 1)
                        #     )
                        #     element_internal_epsi_1[:, epsi_dir] += _w_q_c * (
                        #         element.gradients_operators[_qc].T @ material.mat_data.s1.thermodynamic_forces[_qp]
                        #     )
                        #     mgis_bv.revert(material.mat_data)
                        # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        transformation_gradient = element.compute_transformation_gradient(
                            self.field, _qc, element_unknown_increment
                        )
                        material.mat_data.s1.gradients[_qp] = transformation_gradient
                        # --- INTEGRATE BEHAVIOUR LAW
                        integ_res = mgis_bv.integrate(material.mat_data, material.integration_type, _dt, _qp, (_qp + 1))
                        # --- VOLUMETRIC FORCES
                        v = self.finite_element.cell_basis_l.evaluate_function(
                            _x_q_c, element.cell.shape.centroid, element.cell.shape.diameter
                        )
                        for load in self.loads:
                            vl = _w_q_c * v * load.function(time_step, _x_q_c)
                            _re0 = load.direction * _cl
                            _re1 = (load.direction + 1) * _cl
                            element_external_forces[_re0:_re1] += vl
                        # --- COMPUTE STIFFNESS MATRIX CONTRIBUTION AT QUADRATURE POINT
                        element_stiffness_matrix += _w_q_c * (
                            element.gradients_operators[_qc].T
                            @ material.mat_data.K[_qp]
                            @ element.gradients_operators[_qc]
                        )
                        # --- COMPUTE INTERNAL FORCES CONTRIBUTION AT QUADRATURE POINT
                        element_internal_forces += _w_q_c * (
                            element.gradients_operators[_qc].T @ material.mat_data.s1.thermodynamic_forces[_qp]
                        )
                        _qp += 1
                    # # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    if iteration > 0:
                        numerical_element_stiffness_matrix = (numerical_forces_bs_0[_element_index] - numerical_forces_bs_1[_element_index]) / (
                                2.0 * noise
                        )
                        check_num = np.max(np.abs(numerical_element_stiffness_matrix - element_stiffness_matrix)) / np.max(np.abs(element_stiffness_matrix))
                        if check_num > noise_tolerance:
                            print("AT ELEM : {} | CHECK_BEFORE_STAB : {}".format(str(_element_index).zfill(6), check_num))
                    # # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    # # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    # check_threshold = 10.
                    # numerical_element_stiffness_matrix = (element_internal_epsi_0 - element_internal_epsi_1) / (
                    #     2.0 * noise
                    # )
                    # # check_num = np.max(np.abs(numerical_element_stiffness_matrix)) - np.max(np.abs(element_stiffness_matrix))
                    # check_num = np.max(np.abs(numerical_element_stiffness_matrix - element_stiffness_matrix))
                    # if check_num > check_threshold:
                    #     print("AT ELEM : {} | CHECK_BEFORE_STAB : {}".format(str(_element_index).zfill(6), check_num))
                    # # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    # # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # --- ADDING STABILIZATION CONTRIBUTION AT THE ELEMENT LEVEL
                    element_stiffness_matrix += material.stabilization_parameter * element.stabilization_operator
                    # # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    if iteration > 0:
                        numerical_element_stiffness_matrix = (numerical_forces_as_0[_element_index] - numerical_forces_as_1[_element_index]) / (
                                2.0 * noise
                        )
                        check_num = np.max(np.abs(numerical_element_stiffness_matrix - element_stiffness_matrix)) / np.max(np.abs(element_stiffness_matrix))
                        if check_num > noise_tolerance:
                            print("AT ELEM : {} | CHECK_AFTER_STAB : {}".format(str(_element_index).zfill(6), check_num))
                    # # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    # # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    # for epsi_dir in range(element.element_size):
                    #     element_internal_epsi_0[:, epsi_dir] += (
                    #         material.stabilization_parameter
                    #         * element.stabilization_operator
                    #         @ element_unknown_epsi_0[:, epsi_dir]
                    #     )
                    #     element_internal_epsi_1[:, epsi_dir] += (
                    #         material.stabilization_parameter
                    #         * element.stabilization_operator
                    #         @ element_unknown_epsi_1[:, epsi_dir]
                    #     )
                    # numerical_element_stiffness_matrix = (element_internal_epsi_0 - element_internal_epsi_1) / (
                    #     2.0 * noise
                    # )
                    # # check_num = np.max(np.abs(numerical_element_stiffness_matrix)) - np.max(np.abs(element_stiffness_matrix))
                    # check_num = np.max(np.abs(numerical_element_stiffness_matrix - element_stiffness_matrix))
                    # if check_num > check_threshold:
                    #     print("AT ELEM : {} | CHECK_AFTER_STAB : {}".format(str(_element_index).zfill(6), check_num))
                    # # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    # # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # --- ADDING STABILIZATION CONTRIBUTION AT THE ELEMENT LEVEL
                    element_internal_forces += (
                        material.stabilization_parameter * element.stabilization_operator @ element_unknown_increment
                    )
                    # --- COMPUTE RESIDUAL AFTER VOLUMETRIC CONTRIBUTION
                    element_residual = element_internal_forces - element_external_forces
                    # --- CONDENSATION
                    K_cond, R_cond = element.make_condensation(
                        self.field, self.finite_element, element_stiffness_matrix, element_residual
                    )
                    # --- ASSEMBLY
                    for _i_local, _i_global in enumerate(element.faces_indices):
                        _rg0 = (_dx * _fk) * _i_global
                        _rg1 = (_dx * _fk) * (_i_global + 1)
                        _re0 = _i_local * (_dx * _fk)
                        _re1 = (_i_local + 1) * (_dx * _fk)
                        residual[_rg0:_rg1] += R_cond[_re0:_re1]
                        for _j_local, _j_global in enumerate(element.faces_indices):
                            _cg0 = _j_global * (_dx * _fk)
                            _cg1 = (_j_global + 1) * (_dx * _fk)
                            _ce0 = _j_local * (_dx * _fk)
                            _ce1 = (_j_local + 1) * (_dx * _fk)
                            tangent_matrix[_rg0:_rg1, _cg0:_cg1] += K_cond[_re0:_re1, _ce0:_ce1]
                    # --- SET EXTERNAL FORCES COEFFICIENT
                    if np.max(np.abs(element_external_forces)) > external_forces_coefficient:
                        external_forces_coefficient = np.max(np.abs(element_external_forces))
                # --------------------------------------------------------------------------------------------------
                # FOR FACE LOOP
                # --------------------------------------------------------------------------------------------------
                iter_face_constraint = 0
                for boundary_condition in self.boundary_conditions:
                    # --- DISPLACEMENT CONDITIONS
                    if boundary_condition.boundary_type == BoundaryType.DISPLACEMENT:
                        for element in self.elements:
                            # element_unknown_increment = element.get_element_unknwon_increment(
                            #     self.field, self.finite_element, unknown_increment
                            # )
                            for f_local, f_global in enumerate(element.faces_indices):
                                if (
                                    f_global
                                    in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]
                                ):
                                    _rl0 = _system_size + iter_face_constraint * _fk
                                    _rl1 = _system_size + (iter_face_constraint + 1) * _fk
                                    _re0 = _cl * _dx + (f_local * _dx * _fk) + boundary_condition.direction * _fk
                                    _re1 = _cl * _dx + (f_local * _dx * _fk) + (boundary_condition.direction + 1) * _fk
                                    _rg0 = f_global * _fk * _dx + _fk * boundary_condition.direction
                                    _rg1 = f_global * _fk * _dx + _fk * (boundary_condition.direction + 1)
                                    # -------
                                    # face_displacement = element_unknown_increment[_r0:_r1]
                                    # face_displacement = np.copy(element_unknown_increment[_c0:_c1])
                                    # face_displacement = element_unknown_increment[_c0:_c1]
                                    face_lagrange = unknown_increment[_rl0:_rl1]
                                    face_displacement = unknown_increment[_rg0:_rg1]
                                    _m_psi_psi_face = np.zeros((_fk, _fk), dtype=real)
                                    _v_face_imposed_displacement = np.zeros((_fk,), dtype=real)
                                    for _qf in range(len(element.faces[f_local].quadrature_weights)):
                                        _x_q_f = element.faces[f_local].quadrature_points[:, _qf]
                                        _w_q_f = element.faces[f_local].quadrature_weights[_qf]
                                        v = self.finite_element.face_basis_k.evaluate_function(
                                            _x_q_f,
                                            element.faces[f_local].shape.centroid,
                                            element.faces[f_local].shape.diameter,
                                        )
                                        _v_face_imposed_displacement += (
                                            _w_q_f * v * boundary_condition.function(time_step, _x_q_f)
                                        )
                                        _m_psi_psi_face += blocks.get_face_mass_matrix_in_face(
                                            element.faces[f_local],
                                            self.finite_element.face_basis_k,
                                            self.finite_element.face_basis_k,
                                            _x_q_f,
                                            _w_q_f,
                                        )
                                    _m_psi_psi_face_inv = np.linalg.inv(_m_psi_psi_face)
                                    imposed_face_displacement = _m_psi_psi_face_inv @ _v_face_imposed_displacement
                                    face_displacement_difference = face_displacement - imposed_face_displacement
                                    # --- LAGRANGE INTERNAL FORCES PART
                                    residual[_rg0:_rg1] += normalization_lagrange_coefficient * face_lagrange
                                    # --- LAGRANGE MULTIPLIERS PART
                                    residual[_rl0:_rl1] += (
                                        normalization_lagrange_coefficient * face_displacement_difference
                                    )
                                    # --- LAGRANGE MATRIX PART
                                    tangent_matrix[_rl0:_rl1, _rg0:_rg1] += normalization_lagrange_coefficient * np.eye(
                                        _fk
                                    )
                                    tangent_matrix[_rg0:_rg1, _rl0:_rl1] += normalization_lagrange_coefficient * np.eye(
                                        _fk
                                    )
                                    # --- SET EXTERNAL FORCES COEFFICIENT
                                    lagrange_external_forces = (
                                        normalization_lagrange_coefficient * imposed_face_displacement
                                    )
                                    if np.max(np.abs(lagrange_external_forces)) > external_forces_coefficient:
                                        external_forces_coefficient = np.max(np.abs(lagrange_external_forces))
                                    # if np.max(np.abs(residual[_rg0:_rg1])) > external_forces_coefficient:
                                    #     external_forces_coefficient = np.max(np.abs(residual[_rg0:_rg1]))
                                    # --- ITER ON LAGRANGE CONSTRAINED FACE
                                    iter_face_constraint += 1
                    elif boundary_condition.boundary_type == BoundaryType.PRESSURE:
                        for element in self.elements:
                            # element_unknown_increment = element.get_element_unknwon_increment(
                            #     self.field, self.finite_element, unknown_increment
                            # )
                            for f_local, f_global in enumerate(element.faces_indices):
                                if (
                                    f_global
                                    in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]
                                ):
                                    for qf in range(len(element.faces[f_local].quadrature_weights)):
                                        _x_q_f = element.faces[f_local].quadrature_points[:, qf]
                                        _w_q_f = element.faces[f_local].quadrature_weights[qf]
                                        v = self.finite_element.face_basis_k.evaluate_function(
                                            _x_q_f,
                                            element.faces[f_local].shape.centroid,
                                            element.faces[f_local].shape.diameter,
                                        )
                                        vf = _w_q_f * v * boundary_condition.function(time_step, _x_q_f)
                                        # _c0 = _dx * _cl + f_local * _dx * _fk + boundary_condition.direction * _fk
                                        # _c1 = _dx * _cl + f_local * _dx * _fk + (boundary_condition.direction + 1) * _fk
                                        _rg0 = f_global * _fk * _dx + _fk * boundary_condition.direction
                                        _rg1 = f_global * _fk * _dx + _fk * (boundary_condition.direction + 1)
                                        # residual[_c0:_c1] += vf
                                        residual[_rg0:_rg1] -= vf
                                        # --- SET EXTERNAL FORCES COEFFICIENT
                                        if np.max(np.abs(vf)) > external_forces_coefficient:
                                            external_forces_coefficient = np.max(np.abs(vf))
                # --------------------------------------------------------------------------------------------------
                # RESIDUAL EVALUATION
                # --------------------------------------------------------------------------------------------------
                tol_vect = np.ones((_constrained_system_size,), dtype=real) * self.tolerance
                if external_forces_coefficient == 0.0:
                    external_forces_coefficient = 1.0
                res_eval = (1.0 / external_forces_coefficient) * residual
                print("ITER : {} | RES_MAX : {}".format(str(iteration).zfill(4), max(np.abs(res_eval))))
                residual_values.append(max(np.abs(res_eval)))
                if (np.abs(res_eval) < tol_vect).all():
                    # ----------------------------------------------------------------------------------------------
                    # UPDATE INTERNAL VARIABLES
                    # ----------------------------------------------------------------------------------------------
                    mgis_bv.update(material.mat_data)
                    print("ITERATIONS : {}".format(iteration + 1))
                    self.write_vertex_res_files(file_suffix, unknown_increment)
                    self.write_quadrature_points_res_files(file_suffix, unknown_increment, material)
                    # plt.plot(range(len(residual_values)), residual_values)
                    # plt.show()
                    break
                else:
                    # ----------------------------------------------------------------------------------------------
                    # SOLVE SYSTEM
                    # ----------------------------------------------------------------------------------------------
                    sparse_global_matrix = csr_matrix(-tangent_matrix)
                    correction = spsolve(sparse_global_matrix, residual)
                    unknown_increment += correction
        return

    def get_elements(self, integration_order: int, quadrature_type: QuadratureType = QuadratureType.GAUSS):
        """

        Args:
            integration_order:
            quadrature_type:

        Returns:

        """
        elements = []
        _fk = self.finite_element.face_basis_k.dimension
        _dx = self.field.field_dimension
        _cl = self.finite_element.cell_basis_l.dimension
        for cell_index in range(self.mesh.number_of_cells_in_mesh):
            cell_vertices_connectivity = self.mesh.cells_vertices_connectivity[cell_index]
            cell_faces_connectivity = self.mesh.cells_faces_connectivity[cell_index]
            cell_shape_type = self.mesh.cells_shape_types[cell_index]
            cell_vertices = self.mesh.vertices[:, cell_vertices_connectivity]
            element_cell = Cell(cell_shape_type, cell_vertices, integration_order, quadrature_type=quadrature_type)
            element_faces = []
            element_faces_indices = []
            for global_face_index in cell_faces_connectivity:
                element_faces_indices.append(global_face_index)
                face_vertices_indices = self.mesh.faces_vertices_connectivity[global_face_index]
                face_vertices = self.mesh.vertices[:, face_vertices_indices]
                face_shape_type = self.mesh.faces_shape_types[global_face_index]
                face = Face(face_shape_type, face_vertices, integration_order, quadrature_type=quadrature_type)
                element_faces.append(face)
            element = Element(
                self.field,
                self.finite_element,
                element_cell,
                element_faces,
                element_faces_indices,
                [[(np.nan, np.nan) for j in range(_dx)] for i in range(len(element_faces_indices))],
                [[(np.nan, np.nan) for j in range(_dx)] for i in range(len(element_faces_indices))],
                [[(np.nan, np.nan) for j in range(_dx)] for i in range(len(element_faces_indices))],
            )
            elements.append(element)
            del element_cell
            del element_faces
        constrained_system_size, system_size = systm.get_total_system_size(
            self.field, self.finite_element, self.mesh, self.boundary_conditions
        )
        # iter_face_constraint = 0
        # for element in elements:
        #     for boundary_condition in self.boundary_conditions:
        #         if boundary_condition.boundary_type == BoundaryType.DISPLACEMENT:
        #             for f_local, f_global in enumerate(element.faces_indices):
        #                 if f_global in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]:
        #                     _l0 = system_size + iter_face_constraint * _fk
        #                     _l1 = system_size + (iter_face_constraint + 1) * _fk
        #                     _c0 = _cl * _dx + (f_local * _dx * _fk) + boundary_condition.direction * _fk
        #                     _c1 = _cl * _dx + (f_local * _dx * _fk) + (boundary_condition.direction + 1) * _fk
        #                     _r0 = f_global * _fk * _dx + _fk * boundary_condition.direction
        #                     _r1 = f_global * _fk * _dx + _fk * (boundary_condition.direction + 1)
        #                     element.faces_lagrange_system_row_positions[f_local][boundary_condition.direction] = (_l0, _l1)
        #                     element.faces_lagrange_system_col_positions[f_local][boundary_condition.direction] = (_r0, _r1)
        #                     element.faces_lagrange_local_positions[f_local][boundary_condition.direction] = (_c0, _c1)
        #                     iter_face_constraint += 1
        iter_face_constraint = 0
        for boundary_condition in self.boundary_conditions:
            if boundary_condition.boundary_type == BoundaryType.DISPLACEMENT:
                for element in elements:
                    for f_local, f_global in enumerate(element.faces_indices):
                        if f_global in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]:
                            _l0 = system_size + iter_face_constraint * _fk
                            _l1 = system_size + (iter_face_constraint + 1) * _fk
                            _c0 = _cl * _dx + (f_local * _dx * _fk) + boundary_condition.direction * _fk
                            _c1 = _cl * _dx + (f_local * _dx * _fk) + (boundary_condition.direction + 1) * _fk
                            _r0 = f_global * _fk * _dx + _fk * boundary_condition.direction
                            _r1 = f_global * _fk * _dx + _fk * (boundary_condition.direction + 1)
                            element.faces_lagrange_system_row_positions[f_local][boundary_condition.direction] = (
                                _l0,
                                _l1,
                            )
                            element.faces_lagrange_system_col_positions[f_local][boundary_condition.direction] = (
                                _r0,
                                _r1,
                            )
                            element.faces_lagrange_local_positions[f_local][boundary_condition.direction] = (_c0, _c1)
                            iter_face_constraint += 1
        return elements

    def __check_loads(self, loads: List[Load]):
        """

        Args:
            loads:

        Returns:

        """
        if loads is None:
            return
        if isinstance(loads, list):
            if self.field.field_dimension >= len(loads) > 0:
                for i in range(len(loads)):
                    if isinstance(loads[i], Load):
                        if loads[i].direction < self.field.field_dimension:
                            continue
                        else:
                            raise ValueError
                    else:
                        raise TypeError("loads must be a list of Load")
            else:
                ValueError("loads must be a list of Load of size =< {}".format(self.field.field_dimension))
        else:
            raise TypeError("loads must be a list of Load of size =< {}".format(self.field.field_dimension))
        return

    def __check_boundary_conditions(self, boundary_conditions: List[BoundaryCondition]):
        """

        Args:
            boundary_conditions:

        Returns:

        """
        if isinstance(boundary_conditions, list):
            for boundary_condition in boundary_conditions:
                if isinstance(boundary_condition, BoundaryCondition):
                    if boundary_condition.boundary_name in self.mesh.faces_boundaries_connectivity.keys():
                        if boundary_condition.direction < self.field.field_dimension:
                            continue
                        else:
                            raise ValueError
                    else:
                        raise KeyError
                else:
                    raise TypeError
        else:
            raise TypeError
        return

        # suffix = str(time_step_index).zfill(5)
        # with open(get_res_file_path("vtx_res", suffix), "w") as res_file_g:
        #     for x_dir in range(self.finite_element.euclidean_dimension):
        #         res_file_g.write("X_{} ".format(x_dir))
        #     for x_dir in range(self.finite_element.euclidean_dimension):
        #         res_file_g.write("U_{} ".format(x_dir))
        #     # for vertex_count in self.mesh.get_number_of_vertices_in_mesh()
        #     for vertex_count in range(self.mesh.number_of_vertices_in_mesh):
        #         for x_dir in range(self.finite_element.euclidean_dimension):
        #             v = self.mesh.vertices[x_dir, vertex_count]
        #             res_file_g.write("")

        # def __get_strain_vals(component: int):
        #     strain_vals = [
        #         material.mat_data.s1.gradients[nq, component]
        #         for nq in range(self.mesh.number_of_cell_quadrature_points_in_mesh)
        #     ]
        #     return min(strain_vals), max(strain_vals)
        #
        # def __get_stress_vals(component: int):
        #     stress_vals = [
        #         material.mat_data.s1.thermodynamic_forces[nq, component]
        #         for nq in range(self.mesh.number_of_cell_quadrature_points_in_mesh)
        #     ]
        #     return min(stress_vals), max(stress_vals)
        #
        # def __get_piola_kirchoff_vals(component: int):
        #     # stress_vals = [
        #     #     material.mat_data.s1.thermodynamic_forces[nq, component]
        #     #     for nq in range(self.mesh.number_of_cell_quadrature_points_in_mesh)
        #     # ]
        #     stress_vals = []
        #     for qp in range(material.mat_data.n):
        #         F = np.zeros((3, 3))
        #         F[0, 0] = material.mat_data.s1.gradients[qp, 0]
        #         F[1, 1] = material.mat_data.s1.gradients[qp, 1]
        #         F[2, 2] = material.mat_data.s1.gradients[qp, 2]
        #         F[0, 1] = material.mat_data.s1.gradients[qp, 3]
        #         F[1, 0] = material.mat_data.s1.gradients[qp, 4]
        #         PK = np.zeros((3, 3))
        #         PK[0, 0] = material.mat_data.s1.thermodynamic_forces[qp, 0]
        #         PK[1, 1] = material.mat_data.s1.thermodynamic_forces[qp, 1]
        #         PK[2, 2] = material.mat_data.s1.thermodynamic_forces[qp, 2]
        #         PK[0, 1] = material.mat_data.s1.thermodynamic_forces[qp, 3]
        #         PK[1, 0] = material.mat_data.s1.thermodynamic_forces[qp, 4]
        #         J = np.linalg.det(F)
        #         # F_T_inv = np.linalg.inv(F.T)
        #         sig = (1.0 / J) * PK @ F.T
        #         sig_vect = np.zeros((5,))
        #         sig_vect[0] = sig[0, 0]
        #         sig_vect[1] = sig[1, 1]
        #         sig_vect[2] = sig[2, 2]
        #         sig_vect[3] = sig[0, 1]
        #         sig_vect[4] = sig[1, 0]
        #         stress_vals.append(sig_vect[component])
        #     return min(stress_vals), max(stress_vals)
        #
        # def __write_finite_strain():
        #     res_file.write("{} ".format(time_step))
        #     res_file.write("{} {} ".format(__get_strain_vals(0)[0], __get_strain_vals(0)[1]))
        #     res_file.write("{} {} ".format(__get_strain_vals(1)[0], __get_strain_vals(1)[1]))
        #     res_file.write("{} {} ".format(__get_strain_vals(2)[0], __get_strain_vals(2)[1]))
        #     res_file.write("{} {} ".format(__get_strain_vals(3)[0], __get_strain_vals(3)[1]))
        #     res_file.write("{} {} ".format(__get_strain_vals(4)[0], __get_strain_vals(4)[1]))
        #     # -------------
        #     res_file.write(
        #         "{} {} ".format(__get_piola_kirchoff_vals(0)[0], __get_piola_kirchoff_vals(0)[1])
        #     )
        #     res_file.write(
        #         "{} {} ".format(__get_piola_kirchoff_vals(1)[0], __get_piola_kirchoff_vals(1)[1])
        #     )
        #     res_file.write(
        #         "{} {} ".format(__get_piola_kirchoff_vals(2)[0], __get_piola_kirchoff_vals(2)[1])
        #     )
        #     res_file.write(
        #         "{} {} ".format(__get_piola_kirchoff_vals(3)[0], __get_piola_kirchoff_vals(3)[1])
        #     )
        #     res_file.write(
        #         "{} {} ".format(__get_piola_kirchoff_vals(4)[0], __get_piola_kirchoff_vals(4)[1])
        #     )
        #     res_file.write("\n")
        #
        # def __write_small_strain():
        #     res_file.write("{} ".format(time_step))
        #     res_file.write("{} {} ".format(__get_strain_vals(0)[0], __get_strain_vals(0)[1]))
        #     res_file.write("{} {} ".format(__get_strain_vals(1)[0], __get_strain_vals(1)[1]))
        #     res_file.write("{} {} ".format(__get_strain_vals(2)[0], __get_strain_vals(2)[1]))
        #     res_file.write("{} {} ".format(__get_strain_vals(3)[0], __get_strain_vals(3)[1]))
        #     # -------------
        #     res_file.write("{} {} ".format(__get_stress_vals(0)[0], __get_stress_vals(0)[1]))
        #     res_file.write("{} {} ".format(__get_stress_vals(1)[0], __get_stress_vals(1)[1]))
        #     res_file.write("{} {} ".format(__get_stress_vals(2)[0], __get_stress_vals(2)[1]))
        #     res_file.write("{} {} ".format(__get_stress_vals(3)[0], __get_stress_vals(3)[1]))
        #     res_file.write("\n")

        # __write_finite_strain()
        # __write_small_strain()
        # strain_type = get_strain_type(self.finite_element.field_type)
        # if strain_type == StrainTempType.FINITE_STRAIN:
        #     __write_finite_strain()
        # else:
        #     __write_small_strain()
        # if self.finite_element.field_type in [
        #     FieldType.DISPLACEMENT_FINITE_STRAIN_PLANE_STRESS,
        #     FieldType.DISPLACEMENT_FINITE_STRAIN_PLANE_STRAIN,
        #     FieldType.DISPLACEMENT_FINITE_STRAIN,
        # ]:
        #     __write_finite_strain()
        # else:
        #     __write_small_strain()

        # __plot()


    # def solve_newton_2(self, material: Material, reset_displacement_at_time_step: bool):
    #     """
    #
    #     Args:
    #         material:
    #         reset_displacement_at_time_step:
    #
    #     Returns:
    #
    #     """
    #     self.clean_res_dir()
    #     _dx = self.field.field_dimension
    #     _fk = self.finite_element.face_basis_k.dimension
    #     _cl = self.finite_element.cell_basis_l.dimension
    #     external_forces_coefficient = 1.0
    #     normalization_lagrange_coefficient = material.lagrange_parameter
    #     # ----------------------------------------------------------------------------------------------------------
    #     # SET SYSTEM SIZE
    #     # ----------------------------------------------------------------------------------------------------------
    #     _constrained_system_size, _system_size = systm.get_total_system_size(
    #         self.field, self.finite_element, self.mesh, self.boundary_conditions
    #     )
    #     unknown_increment = np.zeros((_constrained_system_size))
    #     """correction = np.zeros((_constrained_system_size))
    #     """
    #     for time_step_index, time_step in enumerate(self.time_steps):
    #         # --- SET TEMPERATURE
    #         if reset_displacement_at_time_step:
    #             unknown_increment = np.zeros((_constrained_system_size))
    #         material.set_temperature()
    #         # --- PRINT DATA
    #         print("-------------")
    #         print("TIME_STEP : {}".format(time_step))
    #         # --- WRITE RES FILES
    #         residual_values = []
    #         file_suffix = "{}".format(time_step_index).zfill(6)
    #         self.create_vertex_res_files(file_suffix)
    #         self.create_quadrature_points_res_files(file_suffix)
    #         for iteration in range(self.number_of_iterations):
    #             # --------------------------------------------------------------------------------------------------
    #             # SET SYSTEM MATRIX AND VECTOR
    #             # --------------------------------------------------------------------------------------------------
    #             tangent_matrix = np.zeros((_constrained_system_size, _constrained_system_size))
    #             residual = np.zeros((_constrained_system_size))
    #             # --------------------------------------------------------------------------------------------------
    #             # SET TIME INCREMENT
    #             # --------------------------------------------------------------------------------------------------
    #             if time_step_index == 0:
    #                 _dt = time_step
    #             else:
    #                 _dt = time_step - self.time_steps[time_step_index - 1]
    #             _dt = 0.0
    #             # --------------------------------------------------------------------------------------------------
    #             # FOR ELEMENT LOOP
    #             # --------------------------------------------------------------------------------------------------
    #             _qp = 0
    #             iter_face_constraint = 0
    #             for element in self.elements:
    #                 # --- GET ELEMENT UNKNOWN INCREMENT
    #                 element_unknown_increment = element.get_element_unknwon_correction(
    #                     self.field, self.finite_element, unknown_increment
    #                 )
    #                 # --- INITIALIZE MATRIX AND VECTORS
    #                 element_stiffness_matrix = np.zeros((element.element_size, element.element_size))
    #                 element_internal_forces = np.zeros((element.element_size,))
    #                 element_external_forces = np.zeros((element.element_size,))
    #                 # --- RUN OVER EACH QUADRATURE POINT
    #                 for _qc in range(len(element.cell.quadrature_weights)):
    #                     _w_q_c = element.cell.quadrature_weights[_qc]
    #                     _x_q_c = element.cell.quadrature_points[:, _qc]
    #                     # --- COMPUTE STRAINS AND SET THEM IN THE BEHAVIOUR LAW
    #                     transformation_gradient = element.compute_transformation_gradient(
    #                         self.field, _qc, element_unknown_increment
    #                     )
    #                     material.mat_data.s1.gradients[_qp] = transformation_gradient
    #                     # --- INTEGRATE BEHAVIOUR LAW
    #                     integ_res = mgis_bv.integrate(material.mat_data, material.integration_type, _dt, _qp, (_qp + 1))
    #                     # --- VOLUMETRIC FORCES
    #                     v = self.finite_element.cell_basis_l.evaluate_function(
    #                         _x_q_c, element.cell.shape.centroid, element.cell.shape.diameter
    #                     )
    #                     for load in self.loads:
    #                         vl = _w_q_c * v * load.function(time_step, _x_q_c)
    #                         _c0 = load.direction * _cl
    #                         _c1 = (load.direction + 1) * _cl
    #                         element_external_forces[_c0:_c1] += vl
    #                     # --- COMPUTE STIFFNESS MATRIX CONTRIBUTION AT QUADRATURE POINT
    #                     element_stiffness_matrix += _w_q_c * (
    #                         element.gradients_operators[_qc].T
    #                         @ material.mat_data.K[_qp]
    #                         @ element.gradients_operators[_qc]
    #                     )
    #                     # --- COMPUTE STIFFNESS MATRIX CONTRIBUTION AT QUADRATURE POINT
    #                     element_internal_forces += _w_q_c * (
    #                         element.gradients_operators[_qc].T @ material.mat_data.s1.thermodynamic_forces[_qp]
    #                     )
    #                     _qp += 1
    #                 # --- ADDING STABILIZATION CONTRIBUTION AT THE ELEMENT LEVEL
    #                 element_stiffness_matrix += material.stabilization_parameter * element.stabilization_operator
    #                 # --- ADDING STABILIZATION CONTRIBUTION AT THE ELEMENT LEVEL
    #                 element_internal_forces += (
    #                     material.stabilization_parameter * element.stabilization_operator @ element_unknown_increment
    #                 )
    #                 # --- BOUNDARY CONDITIONS
    #                 for boundary_condition in self.boundary_conditions:
    #                     # --- DISPLACEMENT CONDITIONS
    #                     if boundary_condition.boundary_type == BoundaryType.DISPLACEMENT:
    #                         for f_local, f_global in enumerate(element.faces_indices):
    #                             if (
    #                                 f_global
    #                                 in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]
    #                             ):
    #                                 _l0 = _system_size + iter_face_constraint * _fk
    #                                 _l1 = _system_size + (iter_face_constraint + 1) * _fk
    #                                 _c0 = _cl * _dx + (f_local * _dx * _fk) + boundary_condition.direction * _fk
    #                                 _c1 = _cl * _dx + (f_local * _dx * _fk) + (boundary_condition.direction + 1) * _fk
    #                                 _r0 = f_global * _fk * _dx + _fk * boundary_condition.direction
    #                                 _r1 = f_global * _fk * _dx + _fk * (boundary_condition.direction + 1)
    #                                 # -------
    #                                 # face_displacement = element_unknown_increment[_r0:_r1]
    #                                 # face_displacement = np.copy(element_unknown_increment[_c0:_c1])
    #                                 # face_displacement = element_unknown_increment[_c0:_c1]
    #                                 face_lagrange = unknown_increment[_l0:_l1]
    #                                 face_displacement = unknown_increment[_r0:_r1]
    #                                 _m_psi_psi_face = np.zeros((_fk, _fk), dtype=real)
    #                                 _v_face_imposed_displacement = np.zeros((_fk,), dtype=real)
    #                                 for _qf in range(len(element.faces[f_local].quadrature_weights)):
    #                                     _x_q_f = element.faces[f_local].quadrature_points[:, _qf]
    #                                     _w_q_f = element.faces[f_local].quadrature_weights[_qf]
    #                                     v = self.finite_element.face_basis_k.evaluate_function(
    #                                         _x_q_f,
    #                                         element.faces[f_local].shape.centroid,
    #                                         element.faces[f_local].shape.diameter,
    #                                     )
    #                                     _v_face_imposed_displacement += (
    #                                         _w_q_f * v * boundary_condition.function(time_step, _x_q_f)
    #                                     )
    #                                     _m_psi_psi_face += blocks.get_face_mass_matrix_in_face(
    #                                         element.faces[f_local],
    #                                         self.finite_element.face_basis_k,
    #                                         self.finite_element.face_basis_k,
    #                                         _x_q_f,
    #                                         _w_q_f,
    #                                     )
    #                                 _m_psi_psi_face_inv = np.linalg.inv(_m_psi_psi_face)
    #                                 imposed_face_displacement = _m_psi_psi_face_inv @ _v_face_imposed_displacement
    #                                 face_displacement_difference = face_displacement - imposed_face_displacement
    #                                 # --- LAGRANGE INTERNAL FORCES PART
    #                                 element_internal_forces[_c0:_c1] += (
    #                                     normalization_lagrange_coefficient * face_lagrange
    #                                 )
    #                                 # --- LAGRANGE MULTIPLIERS PART
    #                                 residual[_l0:_l1] += (
    #                                     normalization_lagrange_coefficient * face_displacement_difference
    #                                 )
    #                                 # --- LAGRANGE MATRIX PART
    #                                 tangent_matrix[_l0:_l1, _r0:_r1] += normalization_lagrange_coefficient * np.eye(_fk)
    #                                 tangent_matrix[_r0:_r1, _l0:_l1] += normalization_lagrange_coefficient * np.eye(_fk)
    #                                 # --- SET EXTERNAL FORCES COEFFICIENT
    #                                 lagrange_external_forces = (
    #                                     normalization_lagrange_coefficient * imposed_face_displacement
    #                                 )
    #                                 if np.max(np.abs(lagrange_external_forces)) > external_forces_coefficient:
    #                                     external_forces_coefficient = np.max(np.abs(lagrange_external_forces))
    #                                 iter_face_constraint += 1
    #                     elif boundary_condition.boundary_type == BoundaryType.PRESSURE:
    #                         for f_local, f_global in enumerate(element.faces_indices):
    #                             if (
    #                                 f_global
    #                                 in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]
    #                             ):
    #                                 for qf in range(len(element.faces[f_local].quadrature_weights)):
    #                                     _x_q_f = element.faces[f_local].quadrature_points[:, qf]
    #                                     _w_q_f = element.faces[f_local].quadrature_weights[qf]
    #                                     v = self.finite_element.face_basis_k.evaluate_function(
    #                                         _x_q_f,
    #                                         element.faces[f_local].shape.centroid,
    #                                         element.faces[f_local].shape.diameter,
    #                                     )
    #                                     vf = _w_q_f * v * boundary_condition.function(time_step, _x_q_f)
    #                                     _c0 = _dx * _cl + f_local * _dx * _fk + boundary_condition.direction * _fk
    #                                     _c1 = _dx * _cl + f_local * _dx * _fk + (boundary_condition.direction + 1) * _fk
    #                                     _r0 = f_global * _fk * _dx + _fk * boundary_condition.direction
    #                                     _r1 = f_global * _fk * _dx + _fk * (boundary_condition.direction + 1)
    #                                     element_external_forces[_c0:_c1] += vf
    #                 # --- COMPUTE RESIDUAL AFTER VOLUMETRIC CONTRIBUTION
    #                 element_residual = element_internal_forces - element_external_forces
    #                 # --- CONDENSATION
    #                 K_cond, R_cond = element.make_condensation(
    #                     self.field, self.finite_element, element_stiffness_matrix, element_residual
    #                 )
    #                 # --- ASSEMBLY
    #                 for _i_local, _i_global in enumerate(element.faces_indices):
    #                     _rg0 = (_dx * _fk) * _i_global
    #                     _rg1 = (_dx * _fk) * (_i_global + 1)
    #                     _rl0 = _i_local * (_dx * _fk)
    #                     _rl1 = (_i_local + 1) * (_dx * _fk)
    #                     residual[_rg0:_rg1] += R_cond[_rl0:_rl1]
    #                     for _j_local, _j_global in enumerate(element.faces_indices):
    #                         _cg0 = _j_global * (_dx * _fk)
    #                         _cg1 = (_j_global + 1) * (_dx * _fk)
    #                         _cl0 = _j_local * (_dx * _fk)
    #                         _cl1 = (_j_local + 1) * (_dx * _fk)
    #                         tangent_matrix[_rg0:_rg1, _cg0:_cg1] += K_cond[_rl0:_rl1, _cl0:_cl1]
    #                 # --- SET EXTERNAL FORCES COEFFICIENT
    #                 if np.max(np.abs(element_external_forces)) > external_forces_coefficient:
    #                     external_forces_coefficient = np.max(np.abs(element_external_forces))
    #             # --------------------------------------------------------------------------------------------------
    #             # RESIDUAL EVALUATION
    #             # --------------------------------------------------------------------------------------------------
    #             tol_vect = np.ones((_constrained_system_size,), dtype=real) * self.tolerance
    #             if external_forces_coefficient == 0.0:
    #                 external_forces_coefficient = 1.0
    #             res_eval = (1.0 / external_forces_coefficient) * residual
    #             print("ITER : {} | RES_MAX : {}".format(str(iteration).zfill(4), max(np.abs(res_eval))))
    #             residual_values.append(max(np.abs(res_eval)))
    #             if (np.abs(res_eval) < tol_vect).all():
    #                 # ----------------------------------------------------------------------------------------------
    #                 # UPDATE INTERNAL VARIABLES
    #                 # ----------------------------------------------------------------------------------------------
    #                 mgis_bv.update(material.mat_data)
    #                 print("ITERATIONS : {}".format(iteration + 1))
    #                 self.write_vertex_res_files(file_suffix, unknown_increment)
    #                 self.write_quadrature_points_res_files(file_suffix, unknown_increment, material)
    #                 # plt.plot(range(len(residual_values)), residual_values)
    #                 # plt.show()
    #                 break
    #             else:
    #                 # ----------------------------------------------------------------------------------------------
    #                 # SOLVE SYSTEM
    #                 # ----------------------------------------------------------------------------------------------
    #                 sparse_global_matrix = csr_matrix(-tangent_matrix)
    #                 correction = spsolve(sparse_global_matrix, residual)
    #                 unknown_increment += correction
    #     return

    # def solve_newton_0(self, material: Material, reset_displacement_at_time_step: bool):
    #     """
    #
    #     Args:
    #         material:
    #         reset_displacement_at_time_step:
    #
    #     Returns:
    #
    #     """
    #     self.clean_res_dir()
    #     # with open(os.path.join(get_project_path(), "res/res.csv"), "w") as res_file:
    #     _dx = self.field.field_dimension
    #     _fk = self.finite_element.face_basis_k.dimension
    #     _cl = self.finite_element.cell_basis_l.dimension
    #     external_forces_coefficient = 1.0
    #     normalization_lagrange_coefficient = material.lagrange_parameter
    #     # ----------------------------------------------------------------------------------------------------------
    #     # SET SYSTEM SIZE
    #     # ----------------------------------------------------------------------------------------------------------
    #     _constrained_system_size, _system_size = systm.get_total_system_size(
    #         self.field, self.finite_element, self.mesh, self.boundary_conditions
    #     )
    #     unknown_increment = np.zeros((_constrained_system_size))
    #     # correction = np.zeros((_constrained_system_size))
    #     for time_step_index, time_step in enumerate(self.time_steps):
    #         if reset_displacement_at_time_step:
    #             unknown_increment = np.zeros((_constrained_system_size))
    #         material.set_temperature()
    #         print("-------------")
    #         print("TIME_STEP : {}".format(time_step))
    #         residual_values = []
    #         file_suffix = "{}".format(time_step_index).zfill(6)
    #         self.create_vertex_res_files(file_suffix)
    #         self.create_quadrature_points_res_files(file_suffix)
    #         for iteration in range(self.number_of_iterations):
    #             # --------------------------------------------------------------------------------------------------
    #             # SET SYSTEM MATRIX AND VECTOR
    #             # --------------------------------------------------------------------------------------------------
    #             tangent_matrix = np.zeros((_constrained_system_size, _constrained_system_size))
    #             residual = np.zeros((_constrained_system_size))
    #             # --------------------------------------------------------------------------------------------------
    #             # COMPUTE STRAINS
    #             # --------------------------------------------------------------------------------------------------
    #             _qp = 0
    #             for element in self.elements:
    #                 for _qc in range(len(element.cell.quadrature_weights)):
    #                     element_unknown_increment = element.get_element_unknwon_correction(
    #                         self.field,
    #                         self.finite_element,
    #                         unknown_increment
    #                         # self.field, self.finite_element, correction
    #                     )
    #                     transformation_gradient = element.compute_transformation_gradient(
    #                         self.field, _qc, element_unknown_increment
    #                     )
    #                     material.mat_data.s1.gradients[_qp] = transformation_gradient
    #                     # material.mat_data.s1.gradients[_qp] = (
    #                     #     element.gradients_operators[_qc] @ element_unknown_increment
    #                     # )
    #                     _qp += 1
    #             # --------------------------------------------------------------------------------------------------
    #             # INTEGRATE BEHAVIOUR LAW
    #             # --------------------------------------------------------------------------------------------------
    #             if time_step_index == 0:
    #                 _dt = time_step
    #             else:
    #                 _dt = time_step - self.time_steps[time_step_index - 1]
    #             _dt = 0.0
    #             integ_res = mgis_bv.integrate(material.mat_data, material.integration_type, _dt, 0, material.mat_data.n)
    #             # --------------------------------------------------------------------------------------------------
    #             # COMPUTE STIFFNESS MATRIX AND RESIDUAL
    #             # --------------------------------------------------------------------------------------------------
    #             _qp = 0
    #             for element in self.elements:
    #                 # ---------- compute unknown increment
    #                 element_unknown_increment = element.get_element_unknwon_correction(
    #                     self.field,
    #                     self.finite_element,
    #                     unknown_increment
    #                     # self.field, self.finite_element, correction
    #                 )
    #                 # ---------- compute external forces
    #                 element_external_forces = element.get_external_forces(
    #                     self.field,
    #                     self.finite_element,
    #                     time_step,
    #                     self.mesh.faces_boundaries_connectivity,
    #                     boundary_conditions=self.boundary_conditions,
    #                     loads=self.loads,
    #                 )
    #                 # --------- set pbbb forces coefficient
    #                 element_stiffness_matrix = np.zeros((element.element_size, element.element_size))
    #                 element_internal_forces = np.zeros((element.element_size,))
    #                 for _qc in range(len(element.cell.quadrature_weights)):
    #                     _w_q_c = element.cell.quadrature_weights[_qc]
    #                     # ---------- compute stiffness matrix
    #                     element_stiffness_matrix += _w_q_c * (
    #                         element.gradients_operators[_qc].T
    #                         @ material.mat_data.K[_qp]
    #                         @ element.gradients_operators[_qc]
    #                     )
    #                     # ---------- compute internal forces
    #                     element_internal_forces += _w_q_c * (
    #                         element.gradients_operators[_qc].T @ material.mat_data.s1.thermodynamic_forces[_qp]
    #                     )
    #                 _qp += 1
    #                 # ---------- adding the stabilization contribution to the stiffness matrix
    #                 element_stiffness_matrix += material.stabilization_parameter * element.stabilization_operator
    #                 # ---------- adding the stabilization contribution to internal forces
    #                 element_internal_forces += (
    #                     material.stabilization_parameter * element.stabilization_operator @ element_unknown_increment
    #                 )
    #                 # ---------- adding the Lagrange contributions to internal forces
    #                 for _f_local, _f_global in enumerate(element.faces_indices):
    #                     for _direction in range(_dx):
    #                         _pos_global = element.faces_lagrange_system_row_positions[_f_local][_direction]
    #                         _pos_local = element.faces_lagrange_local_positions[_f_local][_direction]
    #                         if not np.isnan(_pos_global[0]):
    #                             _r0 = _pos_global[0]
    #                             _r1 = _pos_global[1]
    #                             _c0 = _pos_local[0]
    #                             _c1 = _pos_local[1]
    #                             lagrange_multiplier_forces = (
    #                                 normalization_lagrange_coefficient
    #                                 * unknown_increment[_r0:_r1]
    #                                 # normalization_lagrange_coefficient * correction[_r0:_r1]
    #                             )
    #                             element_internal_forces[_c0:_c1] += lagrange_multiplier_forces
    #                 # ---------- residual
    #                 element_residual = element_internal_forces - element_external_forces
    #                 # ---------- set external_forces_coefficient
    #                 if max(np.abs(element_external_forces)) > external_forces_coefficient:
    #                     external_forces_coefficient = max(np.abs(element_external_forces))
    #                 # ----------------------------------------------------------------------------------------------
    #                 # CONDENSATION
    #                 # ----------------------------------------------------------------------------------------------
    #                 K_cond, R_cond = element.make_condensation(
    #                     # self.field, self.finite_element, element_stiffness_matrix, element_residual
    #                     self.field,
    #                     self.finite_element,
    #                     element_stiffness_matrix,
    #                     element_residual,
    #                 )
    #                 # ----------------------------------------------------------------------------------------------
    #                 # ASSEMBLY
    #                 # ----------------------------------------------------------------------------------------------
    #                 for _i_local, _i_global in enumerate(element.faces_indices):
    #                     _rg0 = (_dx * _fk) * _i_global
    #                     _rg1 = (_dx * _fk) * (_i_global + 1)
    #                     _rl0 = _i_local * (_dx * _fk)
    #                     _rl1 = (_i_local + 1) * (_dx * _fk)
    #                     residual[_rg0:_rg1] += R_cond[_rl0:_rl1]
    #                     for _j_local, _j_global in enumerate(element.faces_indices):
    #                         _cg0 = _j_global * (_dx * _fk)
    #                         _cg1 = (_j_global + 1) * (_dx * _fk)
    #                         _cl0 = _j_local * (_dx * _fk)
    #                         _cl1 = (_j_local + 1) * (_dx * _fk)
    #                         tangent_matrix[_rg0:_rg1, _cg0:_cg1] += K_cond[_rl0:_rl1, _cl0:_cl1]
    #             # --------------------------------------------------------------------------------------------------
    #             # IMPOSED DISPLACEMENT
    #             # --------------------------------------------------------------------------------------------------
    #             for boundary_condition in self.boundary_conditions:
    #                 if boundary_condition.boundary_type == BoundaryType.DISPLACEMENT:
    #                     for _elem_index, element in enumerate(self.elements):
    #                         for _f_local, _f_global in enumerate(element.faces_indices):
    #                             if (
    #                                 _f_global
    #                                 in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]
    #                             ):
    #                                 _r0 = _cl * _dx + _f_local * _fk * _dx + boundary_condition.direction * _fk
    #                                 _r1 = _cl * _dx + _f_local * _fk * _dx + (boundary_condition.direction + 1) * _fk
    #                                 element_unknown_increment = element.get_element_unknwon_correction(
    #                                     self.field,
    #                                     self.finite_element,
    #                                     unknown_increment
    #                                     # self.field, self.finite_element, correction
    #                                 )
    #                                 current_face_component_displacement = element_unknown_increment[_r0:_r1]
    #                                 _m_psi_psi_face = np.zeros((_fk, _fk), dtype=real)
    #                                 _v_face_imposed_displacement = np.zeros((_fk,), dtype=real)
    #                                 for _qf in range(len(element.faces[_f_local].quadrature_weights)):
    #                                     _x_q_f = element.faces[_f_local].quadrature_points[:, _qf]
    #                                     _w_q_f = element.faces[_f_local].quadrature_weights[_qf]
    #                                     v = self.finite_element.face_basis_k.evaluate_function(
    #                                         _x_q_f,
    #                                         element.faces[_f_local].shape.centroid,
    #                                         element.faces[_f_local].shape.diameter,
    #                                     )
    #                                     _v_face_imposed_displacement += (
    #                                         _w_q_f * v * boundary_condition.function(time_step, _x_q_f)
    #                                     )
    #                                     _m_psi_psi_face += blocks.get_face_mass_matrix_in_face(
    #                                         element.faces[_f_local],
    #                                         self.finite_element.face_basis_k,
    #                                         self.finite_element.face_basis_k,
    #                                         _x_q_f,
    #                                         _w_q_f,
    #                                     )
    #                                 _m_psi_psi_face_inv = np.linalg.inv(_m_psi_psi_face)
    #                                 imposed_displacement_component_projection = (
    #                                     _m_psi_psi_face_inv @ _v_face_imposed_displacement
    #                                 )
    #                                 _r0 = element.faces_lagrange_system_row_positions[_f_local][
    #                                     boundary_condition.direction
    #                                 ][0]
    #                                 _r1 = element.faces_lagrange_system_row_positions[_f_local][
    #                                     boundary_condition.direction
    #                                 ][1]
    #                                 _c0 = _f_global * _dx * _fk + boundary_condition.direction * _fk
    #                                 _c1 = (_f_global * _dx * _fk) + (boundary_condition.direction + 1) * _fk
    #                                 _local_lagrange_matrix = normalization_lagrange_coefficient * np.eye(_fk)
    #                                 tangent_matrix[_r0:_r1, _c0:_c1] += _local_lagrange_matrix
    #                                 tangent_matrix[_c0:_c1, _r0:_r1] += _local_lagrange_matrix
    #                                 residual[_r0:_r1] += normalization_lagrange_coefficient * (
    #                                     current_face_component_displacement - imposed_displacement_component_projection
    #                                 )
    #                                 if (
    #                                     max(np.abs(imposed_displacement_component_projection))
    #                                     * normalization_lagrange_coefficient
    #                                     > external_forces_coefficient
    #                                 ):
    #                                     external_forces_coefficient = max(
    #                                         np.abs(imposed_displacement_component_projection)
    #                                         * normalization_lagrange_coefficient
    #                                     )
    #             # --------------------------------------------------------------------------------------------------
    #             # RESIDUAL EVALUATION
    #             # --------------------------------------------------------------------------------------------------
    #             tol_vect = np.ones((_constrained_system_size,), dtype=real) * self.tolerance
    #             if external_forces_coefficient == 0.0:
    #                 external_forces_coefficient = 1.0
    #             res_eval = (1.0 / external_forces_coefficient) * residual
    #             # resiudal_noise = np.random.random_sample() * 1.e-6 * max(unknown_increment[:_system_size])
    #             # resiudal_noise = 1.e-2 * max(unknown_increment[:_system_size])
    #             # resiudal_noise_vector = np.full((_constrained_system_size,), resiudal_noise)
    #             # noisy_residual = (residual + resiudal_noise_vector) - (residual - resiudal_noise_vector) / 2.*resiudal_noise
    #             # dirrrection = 0
    #             # for f_index in range(self.mesh.number_of_faces_in_mesh):
    #             #     _c000 = f_index * _fk * _dx + _fk * dirrrection
    #             #     _c001 = f_index * _fk * _dx + _fk * (dirrrection+1)
    #             #     residual[_c000:_c001] = noisy_residual[_c000:_c001]
    #             # residual = noisy_residual
    #             # print("ITER : {} | RES_MAX : {} | NOISE : {}".format(str(iteration).zfill(4), max(np.abs(res_eval)), resiudal_noise))
    #             print("ITER : {} | RES_MAX : {}".format(str(iteration).zfill(4), max(np.abs(res_eval))))
    #             residual_values.append(max(np.abs(res_eval)))
    #             if (np.abs(res_eval) < tol_vect).all():
    #                 # ----------------------------------------------------------------------------------------------
    #                 # UPDATE INTERNAL VARIABLES
    #                 # ----------------------------------------------------------------------------------------------
    #                 mgis_bv.update(material.mat_data)
    #                 print("ITERATIONS : {}".format(iteration + 1))
    #                 self.write_vertex_res_files(file_suffix, unknown_increment)
    #                 self.write_quadrature_points_res_files(file_suffix, unknown_increment, material)
    #                 # plt.plot(range(len(residual_values)), residual_values)
    #                 # plt.show()
    #                 break
    #             else:
    #                 # ----------------------------------------------------------------------------------------------
    #                 # SOLVE SYSTEM
    #                 # ----------------------------------------------------------------------------------------------
    #                 sparse_global_matrix = csr_matrix(-tangent_matrix)
    #                 correction = spsolve(sparse_global_matrix, residual)
    #                 unknown_increment += correction
    #     return

    # def solve_newton(self, material: Material):
    #     for time_step_index, time_step in enumerate(self.time_steps):
    #         material.set_temperature()
    #         print("-------------")
    #         print("TIME_STEP : {}".format(time_step))
    #         for iteration in range(self.iterations):
    #             qp = 0
    #             for element in self.elements:
    #                 for qc in range(len(element.cell.quadrature_weights)):
    #                     elem_incr_unknown = element.get_element_unknwon_increment(
    #                         self.finite_element, self.system.displacement_increment
    #                     )
    #                     material.mat_data.s1.gradients[qp] = element.gradients_operators[qc] @ elem_incr_unknown
    #                     qp += 1
    #             if time_step_index == 0:
    #                 dt = time_step
    #             else:
    #                 dt = time_step - self.time_steps[time_step_index - 1]
    #             dt = 0.0
    #             mgis_bv.integrate(material.mat_data, material.integration_type, dt, 0, material.mat_data.n)
    #             qp = 0
    #             coef_try = 0.0
    #             for element in self.elements:
    #                 element_external_forces = element.get_external_forces(
    #                     self.finite_element,
    #                     time_step,
    #                     self.mesh.faces_boundaries_connectivity,
    #                     boundary_conditions=self.boundary_conditions,
    #                     loads=self.loads,
    #                 )
    #                 element_internal_forces = np.zeros((element.element_size,))
    #                 element_stiffness_matrix = np.zeros((element.element_size, element.element_size))
    #                 for qc in range(len(element.cell.quadrature_weights)):
    #                     w_q_c = element.cell.quadrature_weights[qc]
    #                     # ---------- stiffness matrix
    #                     element_stiffness_matrix += w_q_c * (
    #                         element.gradients_operators[qc].T
    #                         @ material.mat_data.K[qp]
    #                         @ element.gradients_operators[qc]
    #                     )
    #                     # ---------- internal forces
    #                     element_internal_forces += (
    #                         element.gradients_operators[qc].T @ material.mat_data.s1.thermodynamic_forces[qp]
    #                     )
    #                     qp += 1
    #                 element_stiffness_matrix += material.stabilization_parameter * element.stabilization_operator
    #                 element_residual = element_internal_forces - element_external_forces
    #                 # K_cond, Fe_cond = element.make_condensation(
    #                 #     self.finite_element, element_stiffness_matrix, element_external_forces
    #                 # )
    #                 # K_cond, Fe_cond, Fi_cond = element.make_condensation(
    #                 #     self.finite_element, element_stiffness_matrix, element_external_forces, element_internal_forces
    #                 # )
    #                 K_cond, R_cond = element.make_condensation(
    #                     self.finite_element, element_stiffness_matrix, element_residual
    #                 )
    #                 _, Fe_cond = element.make_condensation(
    #                     self.finite_element, element_stiffness_matrix, element_external_forces
    #                 )
    #                 _, Fi_cond = element.make_condensation(
    #                     self.finite_element, element_stiffness_matrix, element_internal_forces
    #                 )
    #                 # self.assemble_element(element, K_cond, Fe_cond)
    #                 # --------------------------------------------------------------------------------------------------
    #                 # ASSEMBLY
    #                 # --------------------------------------------------------------------------------------------------
    #                 fk = self.finite_element.face_basis_k.dimension
    #                 dx = self.finite_element.field_dimension
    #                 for i_local, i_global in enumerate(element.faces_indices):
    #                     rg0 = (dx * fk) * i_global
    #                     rg1 = (dx * fk) * (i_global + 1)
    #                     rl0 = i_local * (dx * fk)
    #                     rl1 = (i_local + 1) * (dx * fk)
    #                     self.system.residual[rg0:rg1] += R_cond[rl0:rl1]
    #                     self.system.internal_forces[rg0:rg1] += Fi_cond[rl0:rl1]
    #                     self.system.external_forces[rg0:rg1] += Fe_cond[rl0:rl1]
    #                     for j_local, j_global in enumerate(element.faces_indices):
    #                         cg0 = j_global * (dx * fk)
    #                         cg1 = (j_global + 1) * (dx * fk)
    #                         cl0 = j_local * (dx * fk)
    #                         cl1 = (j_local + 1) * (dx * fk)
    #                         self.system.tangent_matrix[rg0:rg1, cg0:cg1] += K_cond[rl0:rl1, cl0:cl1]
    #                 # self.assemble_element(element, K_cond, R_cond)
    #             self.system.residual2 = self.system.internal_forces - self.system.external_forces
    #             self.apply_displacement_conditions(time_step, self.boundary_conditions)
    #             # -----
    #             # self.system.internal_forces = self.system.global_matrix @ self.system.displacement_increment
    #             # res = self.system.internal_forces - self.system.global_vector
    #             tol_vect = np.ones((self.system.constrained_system_size,), dtype=real) * self.tolerance
    #             coef = max(np.abs(self.system.external_forces))
    #             if coef == 0.0:
    #                 coef = 1.0
    #             res_eval = (1.0 / coef) * self.system.residual
    #             if (np.abs(res_eval) < tol_vect).all():
    #                 self.system.displacement = self.system.displacement_increment
    #                 mgis_bv.update(material.mat_data)
    #
    #                 strain_vals = [
    #                     material.mat_data.s1.gradients[nq, 0]
    #                     for nq in range(self.mesh.number_of_cell_quadrature_points_in_mesh)
    #                 ]
    #                 min_eps = min(strain_vals)
    #                 max_eps = max(strain_vals)
    #                 print("EPS_MIN : {}".format(min_eps))
    #                 print("EPS_MAX : {}".format(max_eps))
    #                 stress_vals = [
    #                     material.mat_data.s1.thermodynamic_forces[nq, 0]
    #                     for nq in range(self.mesh.number_of_cell_quadrature_points_in_mesh)
    #                 ]
    #                 min_sig = min(stress_vals)
    #                 max_sig = max(stress_vals)
    #                 print("SIG_MIN : {}".format(min_sig))
    #                 print("SIG_MAX : {}".format(max_sig))
    #
    #                 # __plot()
    #                 break
    #             else:
    #                 # self.system.residual = res
    #                 self.system.solve()
    #                 self.system.reset_system()
    #                 self.system.displacement_increment += self.system.correction
    #     return

    # def __plot():
    #
    #     fig, ax0d = plt.subplots(nrows=1, ncols=1)
    #     X, Y = self.mesh.vertices
    #     quadrature_points = []
    #     # div_at_quadrature_points = []
    #     data = np.zeros(self.mesh.number_of_cell_quadrature_points_in_mesh)
    #     qp_count = 0
    #     for element in self.elements:
    #         for qc in range(len(element.cell.quadrature_weights)):
    #             qp = element.cell.quadrature_points[:, qc]
    #             quadrature_points.append(qp)
    #             # div_at_quadrature_point = (
    #             #     element.material_points.strains[qc, 0] + element.material_points.strains[qc, 1]
    #             # )
    #             # div_at_quadrature_points.append(div_at_quadrature_point)
    #             # data[qp_count] = material.mat_data.s1.thermodynamic_forces[qp_count, 1]
    #             data[qp_count] = (
    #                     material.mat_data.s1.gradients[qp_count, 0]
    #                     + material.mat_data.s1.gradients[qp_count, 1]
    #             )
    #             # data[qp_count] = material.mat_data.s1.thermodynamic_forces[qp_count, 0]
    #             # data[qp_count] = material.mat_data.s1.gradients[qp_count, 0]
    #             qp_count += 1
    #     quadrature_points = np.array(quadrature_points).T
    #     x, y = quadrature_points
    #     colors = [(0, 0, 1), (0, 1, 1), (0, 1, 0), (1, 1, 0), (1, 0, 0)]
    #     perso = LinearSegmentedColormap.from_list("perso", colors, N=1000)
    #     vmin = min(data[:])
    #     vmax = max(data[:])
    #     levels = np.linspace(vmin, vmax, 1000, endpoint=True)
    #     ticks = np.linspace(vmin, vmax, 10, endpoint=True)
    #     datad = ax0d.tricontourf(x, y, data[:], cmap=perso, levels=levels)
    #     ax0d.get_xaxis().set_visible(False)
    #     ax0d.get_yaxis().set_visible(False)
    #     ax0d.set_xlabel("map of the domain $\Omega$")
    #     cbar = fig.colorbar(datad, ax=ax0d, ticks=ticks)
    #     cbar.set_label("trace eps", rotation=270, labelpad=15.0)
    #     # plt.savefig("/home/dsiedel/Projects/pythhon/plots/{}.png".format(time_step))
    #     plt.show()
    #     return

    # def assemble_element(self, element: Element, matrix_cond: ndarray, vector_cond: ndarray):
    #     fk = self.finite_element.face_basis_k.dimension
    #     dx = self.finite_element.field_dimension
    #     for i_local, i_global in enumerate(element.faces_indices):
    #         rg0 = (dx * fk) * i_global
    #         rg1 = (dx * fk) * (i_global + 1)
    #         rl0 = i_local * (dx * fk)
    #         rl1 = (i_local + 1) * (dx * fk)
    #         self.system.residual[rg0:rg1] += vector_cond[rl0:rl1]
    #         for j_local, j_global in enumerate(element.faces_indices):
    #             cg0 = j_global * (dx * fk)
    #             cg1 = (j_global + 1) * (dx * fk)
    #             cl0 = j_local * (dx * fk)
    #             cl1 = (j_local + 1) * (dx * fk)
    #             self.system.tangent_matrix[rg0:rg1, cg0:cg1] += matrix_cond[rl0:rl1, cl0:cl1]
    #     return

    # def apply_displacement_conditions(self, time_step: float, boundary_conditions: List[BoundaryCondition] = None):
    #     fe = self.finite_element
    #     dx = fe.field_dimension
    #     fk = fe.face_basis_k.dimension
    #     iter_face_constraint = 0
    #     iter_cons_constraint = 0
    #     norm_ceof = np.max(self.system.tangent_matrix)
    #     for boundary_condition in boundary_conditions:
    #         if boundary_condition.boundary_type == BoundaryType.DISPLACEMENT:
    #             for element in self.elements:
    #                 for f_elemnt, f_global in enumerate(element.faces_indices):
    #                     if f_global in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]:
    #                         m_psi_psi_face = np.zeros((fk, fk), dtype=real)
    #                         vf = np.zeros((fk,), dtype=real)
    #                         for qf in range(len(element.faces[f_elemnt].quadrature_weights)):
    #                             x_q_f = element.faces[f_elemnt].quadrature_points[:, qf]
    #                             w_q_f = element.faces[f_elemnt].quadrature_weights[qf]
    #                             v = fe.face_basis_k.evaluate_function(
    #                                 x_q_f,
    #                                 element.faces[f_elemnt].shape.centroid,
    #                                 element.faces[f_elemnt].shape.diameter,
    #                             )
    #                             vf += w_q_f * v * boundary_condition.function(time_step, x_q_f)
    #                             m_psi_psi_face += blocks.get_face_mass_matrix_in_face(
    #                                 element.faces[f_elemnt], fe.face_basis_k, fe.face_basis_k, x_q_f, w_q_f
    #                             )
    #                         m_psi_psi_face_inv = np.linalg.inv(m_psi_psi_face)
    #                         vo = norm_ceof * m_psi_psi_face_inv @ vf
    #                         l0 = self.system.system_size + iter_face_constraint * fk
    #                         l1 = self.system.system_size + (iter_face_constraint + 1) * fk
    #                         c0 = f_global * dx * fk + boundary_condition.direction * fk
    #                         c1 = (f_global * dx * fk) + (boundary_condition.direction + 1) * fk
    #                         self.system.tangent_matrix[l0:l1, c0:c1] += norm_ceof * np.eye(fk)
    #                         self.system.tangent_matrix[c0:c1, l0:l1] += self.system.tangent_matrix[l0:l1, c0:c1].T
    #                         self.system.residual[l0:l1] += vo
    #                         self.system.residual2[l0:l1] += vo
    #                         self.system.external_forces[l0:l1] += vo
    #                         iter_face_constraint += 1
    #         elif boundary_condition.boundary_type == BoundaryType.SLIDE:
    #             for element in self.elements:
    #                 for f_elemnt, f_global in enumerate(element.faces_indices):
    #                     if f_global in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]:
    #                         vf = boundary_condition.function(time_step, element.faces[0].quadrature_points[:, 0])
    #                         l0 = (
    #                             self.system.system_size
    #                             + self.system.constrained_faces_matrix_size
    #                             + iter_cons_constraint
    #                         )
    #                         c0 = f_global * dx * fk + boundary_condition.direction * fk
    #                         self.system.tangent_matrix[l0, c0] += 1.0
    #                         self.system.tangent_matrix[c0, l0] += 1.0
    #                         self.system.residual[l0] += vf
    #                         # self.system.residual2[l0:l1] += vo
    #                         # self.system.external_forces[l0:l1] += vo
    #                         iter_cons_constraint += 1
    #     return

    # def __plot():
    #
    #     fig, ax0d = plt.subplots(nrows=1, ncols=1)
    #     X, Y = self.mesh.vertices
    #     quadrature_points = []
    #     # div_at_quadrature_points = []
    #     data = np.zeros(self.mesh.number_of_cell_quadrature_points_in_mesh)
    #     qp_count = 0
    #     for element in self.elements:
    #         for qc in range(len(element.cell.quadrature_weights)):
    #             qp = element.cell.quadrature_points[:, qc]
    #             quadrature_points.append(qp)
    #             # div_at_quadrature_point = (
    #             #     element.material_points.strains[qc, 0] + element.material_points.strains[qc, 1]
    #             # )
    #             # div_at_quadrature_points.append(div_at_quadrature_point)
    #             # data[qp_count] = material.mat_data.s1.thermodynamic_forces[qp_count, 1]
    #             data[qp_count] = (
    #                     material.mat_data.s1.gradients[qp_count, 0]
    #                     + material.mat_data.s1.gradients[qp_count, 1]
    #             )
    #             # data[qp_count] = material.mat_data.s1.thermodynamic_forces[qp_count, 0]
    #             data[qp_count] = material.mat_data.s1.gradients[qp_count, 0]
    #             qp_count += 1
    #     quadrature_points = np.array(quadrature_points).T
    #     x, y = quadrature_points
    #     colors = [(0, 0, 1), (0, 1, 1), (0, 1, 0), (1, 1, 0), (1, 0, 0)]
    #     perso = LinearSegmentedColormap.from_list("perso", colors, N=1000)
    #     vmin = min(data[:])
    #     vmax = max(data[:])
    #     levels = np.linspace(vmin, vmax, 1000, endpoint=True)
    #     ticks = np.linspace(vmin, vmax, 10, endpoint=True)
    #     datad = ax0d.tricontourf(x, y, data[:], cmap=perso, levels=levels)
    #     ax0d.get_xaxis().set_visible(False)
    #     ax0d.get_yaxis().set_visible(False)
    #     ax0d.set_xlabel("map of the domain $\Omega$")
    #     cbar = fig.colorbar(datad, ax=ax0d, ticks=ticks)
    #     cbar.set_label("stress yy", rotation=270, labelpad=15.0)
    #     # plt.savefig("/home/dsiedel/Projects/pythhon/plots/{}.png".format(time_step))
    #     plt.show()
    #     return

    # def assemble_system(self, time_step):
    #     elem_count = 0
    #     for element in self.elements:
    #         external_forces = element.get_external_forces(
    #             self.finite_element,
    #             time_step,
    #             self.mesh.faces_boundaries_connectivity,
    #             boundary_conditions=self.boundary_conditions,
    #             loads=self.loads,
    #         )
    #         element_stiffness_matrix = element.get_stiffness_matrix()
    #         K_cond, Fe_cond = element.make_condensation(self.finite_element, external_forces, element_stiffness_matrix)
    #         element.assemble_element2(self.finite_element, Fe_cond, K_cond, self.system)
    #         # ----------------------------------------------------------------------
    #         # internal_forces = element.get_internal_forces2(self.finite_element, self.system.increment)
    #         # external_forces = element.get_external_forces(
    #         #     self.finite_element,
    #         #     time_step,
    #         #     self.mesh.faces_boundaries_connectivity,
    #         #     boundary_conditions=self.boundary_conditions,
    #         #     loads=self.loads,
    #         # )
    #         # element_stiffness_matrix = element.get_stiffness_matrix()
    #         # element_residual = internal_forces - external_forces
    #         # element.assemble_element(self.finite_element, element_residual, element_stiffness_matrix, self.system)
    #         # elem_count += 1
    #     return
    #
    # def assign_strain(self):
    #     for element in self.elements:
    #         # element.set_strain(self.finite_element, self.system.increment2)
    #         s = self.system.system_size
    #         element.set_strain(self.finite_element, self.system.displacement[:s])
    #     return

    # def solve_newton(self,
    #                  material: Material
    #                  ):
    #     # behaviour = mgis_bv.load(
    #     #     material.lib_path, material.behaviour_name, material.hypothesis
    #     # )
    #     # material data manager
    #     # nq = self.mesh.number_of_cell_quadrature_points_in_mesh
    #     # mat_data = mgis_bv.MaterialDataManager(behaviour, nq)
    #     for time_step in self.time_steps:
    #         # m.s0.setMaterialProperty("YoungModulus", 1.0)
    #         # for key, val in material.material_params.items():
    #         #     mat_data.s0.setMaterialProperty(key, val)
    #         #     mat_data.s1.setMaterialProperty(key, val)
    #         # for key, val in material.external_params.items():
    #         #     T = val * np.ones(nq)
    #         #     mgis_bv.setExternalStateVariable(mat_data.s0, key, T, material.external_variable_storage_mode)
    #         #     mgis_bv.setExternalStateVariable(mat_data.s1, key, T, material.external_variable_storage_mode)
    #         # T = temperature * np.ones(nq)
    #         # mgis_bv.setExternalStateVariable(mat_data.s0, "Temperature", T, mgis_bv.MaterialStateManagerStorageMode.ExternalStorage)
    #         # mgis_bv.setExternalStateVariable(mat_data.s1, "Temperature", T, mgis_bv.MaterialStateManagerStorageMode.ExternalStorage)
    #         material.set_temperature()
    #         for iteration in range(self.iterations):
    #             # print("--------")
    #             # print("TIME STEP : {}".format(time_step))
    #             # print("ITERATION : {}".format(iteration))
    #             qp_count = 0
    #             for element in self.elements:
    #                 for qp in range(len(element.cell.quadrature_weights)):
    #                     s = self.system.system_size
    #                     elem_incr = element.get_element_increment(
    #                         self.finite_element, self.system.increment[:s]
    #                     )
    #                     material.mat_data.s0.gradients[qp_count] = (
    #                         element.gradients_operators[qp] @ elem_incr
    #                     )
    #                     qp_count += 1
    #             dt = 0.0
    #             mgis_bv.integrate(material.mat_data, material.integration_type, dt, 0, material.mat_data.n)
    #             qp_count = 0
    #             for element in self.elements:
    #                 external_forces = element.get_external_forces(
    #                     self.finite_element,
    #                     time_step,
    #                     self.mesh.faces_boundaries_connectivity,
    #                     boundary_conditions=self.boundary_conditions,
    #                     loads=self.loads,
    #                 )
    #                 # element_stiffness_matrix = element.get_stiffness_matrix()
    #                 element_stiffness_matrix = np.zeros(
    #                     (element.element_size, element.element_size)
    #                 )
    #                 for qc in range(len(element.cell.quadrature_weights)):
    #                     w_q_c = element.cell.quadrature_weights[qc]
    #                     element_stiffness_matrix += w_q_c * (
    #                             element.gradients_operators[qc].T
    #                             @ material.mat_data.K[qp_count]
    #                             @ element.gradients_operators[qc]
    #                     )
    #                     qp_count += 1
    #                 element_stiffness_matrix += element.stabilization_operator
    #                 K_cond, Fe_cond = element.make_condensation(
    #                     self.finite_element, external_forces, element_stiffness_matrix
    #                 )
    #                 element.assemble_element2(
    #                     self.finite_element, Fe_cond, K_cond, self.system
    #                 )
    #             # self.assemble_system(time_step)
    #             self.apply_displacement_conditions(time_step, self.boundary_conditions)
    #             # -----
    #             s = self.system.system_size
    #             self.system.internal_forces = (
    #                 self.system.global_matrix @ self.system.increment
    #             )[:s]
    #             res = self.system.internal_forces - self.system.external_forces
    #             tol_vect = np.ones((s,), dtype=real) * self.tolerance
    #             if iteration == 0:
    #                 coef = max(np.abs(self.system.external_forces))
    #                 if coef == 0.0:
    #                     coef = 1.0
    #             print("COEF : {}".format(coef))
    #
    #             def __plot():
    #
    #                 fig, ax0d = plt.subplots(nrows=1, ncols=1)
    #                 X, Y = self.mesh.vertices
    #                 quadrature_points = []
    #                 div_at_quadrature_points = []
    #                 for element in self.elements:
    #                     for qc in range(len(element.cell.quadrature_weights)):
    #                         qp = element.cell.quadrature_points[:, qc]
    #                         quadrature_points.append(qp)
    #                         div_at_quadrature_point = (
    #                             element.material_points.strains[qc, 0]
    #                             + element.material_points.strains[qc, 1]
    #                         )
    #                         div_at_quadrature_points.append(div_at_quadrature_point)
    #                 quadrature_points = np.array(quadrature_points).T
    #                 div_at_quadrature_points = np.array(div_at_quadrature_points)
    #                 x, y = quadrature_points
    #                 colors = [(0, 0, 1), (0, 1, 1), (0, 1, 0), (1, 1, 0), (1, 0, 0)]
    #                 perso = LinearSegmentedColormap.from_list("perso", colors, N=1000)
    #                 vmin = min(div_at_quadrature_points[:])
    #                 vmax = max(div_at_quadrature_points[:])
    #                 levels = np.linspace(vmin, vmax, 1000, endpoint=True)
    #                 ticks = np.linspace(vmin, vmax, 10, endpoint=True)
    #                 datad = ax0d.tricontourf(
    #                     x, y, div_at_quadrature_points[:], cmap=perso, levels=levels
    #                 )
    #                 ax0d.get_xaxis().set_visible(False)
    #                 ax0d.get_yaxis().set_visible(False)
    #                 ax0d.set_xlabel("map of the domain $\Omega$")
    #                 cbar = fig.colorbar(datad, ax=ax0d, ticks=ticks)
    #                 cbar.set_label(
    #                     "divergence of the displacement", rotation=270, labelpad=15.0
    #                 )
    #                 # plt.savefig("/home/dsiedel/Projects/pythhon/plots/{}.png".format(time_step))
    #                 plt.show()
    #                 return
    #
    #             if (np.abs(res) < coef * tol_vect).all():
    #                 # if ((1.0 / coef) * np.abs(self.system.global_vector[: self.system.system_size]) < tol_vect).all():
    #                 # ---
    #                 # self.system.displacement += self.system.increment
    #                 # ---
    #                 self.system.displacement = self.system.increment
    #                 # ---
    #                 self.assign_strain()
    #                 # ---
    #                 # res = np.zeros((s,), dtype=real)
    #                 # ---
    #                 # self.system.internal_forces = np.zeros((s,), dtype=real)
    #                 # ---
    #                 # self.system.increment = self.system.displacement
    #                 # ---
    #                 # self.system.increment = np.zeros((self.system.constrained_system_size,), dtype=real)
    #                 mgis_bv.update(material.mat_data)
    #
    #                 __plot()
    #                 # mgis update
    #                 break
    #             else:
    #                 self.system.global_vector[:s] = res
    #                 self.system.solve()
    #                 self.system.reset_system()
    #                 self.system.increment += self.system.correction
    #                 # self.system.increment = self.system.correction
    #                 # ----------------
    #                 # self.system.internal_forces = (self.system.global_matrix @ self.system.correction)[:s]
    #                 # self.system.internal_forces = (self.system.global_matrix @ self.system.increment)[:s]
    #                 # self.assign_strain()
    #                 # self.system.reset_system()
    #                 # mgis integrate
    #     return
    #
    # # def get_m_things(self):
    # #     m_list = []
    # #     for element in self.elements:
    # #         lib = "/home/dsiedel/Projects/pythhon/behaviour/src/libBehaviour.so"
    # #         h = mgis_bv.Hypothesis.PLANESTRAIN
    # #         b = mgis_bv.load(lib, "Elasticity", h)
    # #         nq = len(element.cell.quadrature_weights)
    # #         m = mgis_bv.MaterialDataManager(b, nq)
    # #         m.s0.setMaterialProperty("YoungModulus", 1.124999981250001)
    # #         m.s0.setMaterialProperty("PoissonRatio", 0.499999999999999)
    # #         m.s1.setMaterialProperty("YoungModulus", 1.124999981250001)
    # #         m.s1.setMaterialProperty("PoissonRatio", 0.499999999999999)
    # #         T = 293.15 * np.ones(nq)
    # #         Ts = mgis_bv.MaterialStateManagerStorageMode.ExternalStorage
    # #         mgis_bv.setExternalStateVariable(m.s0, "Temperature", T, Ts)
    # #         mgis_bv.setExternalStateVariable(m.s1, "Temperature", T, Ts)
    # #         m_list.append(m)
    # #     return m_list

    # def __do():
    #     strain_vals = [
    #         material.mat_data.s1.gradients[nq, 0]
    #         for nq in range(self.mesh.number_of_cell_quadrature_points_in_mesh)
    #     ]
    #     min_eps = min(strain_vals)
    #     max_eps = max(strain_vals)
    #     print("EPS_XX_MIN : {}".format(min_eps))
    #     print("EPS_XX_MAX : {}".format(max_eps))
    #     strain_vals = [
    #         material.mat_data.s1.gradients[nq, 1]
    #         for nq in range(self.mesh.number_of_cell_quadrature_points_in_mesh)
    #     ]
    #     min_eps = min(strain_vals)
    #     max_eps = max(strain_vals)
    #     print("EPS_YY_MIN : {}".format(min_eps))
    #     print("EPS_YY_MAX : {}".format(max_eps))
    #     strain_vals = [
    #         material.mat_data.s1.gradients[nq, 2]
    #         for nq in range(self.mesh.number_of_cell_quadrature_points_in_mesh)
    #     ]
    #     min_eps = min(strain_vals)
    #     max_eps = max(strain_vals)
    #     print("EPS_ZZ_MIN : {}".format(min_eps))
    #     print("EPS_ZZ_MAX : {}".format(max_eps))
    #     strain_vals = [
    #         material.mat_data.s1.gradients[nq, 3]
    #         for nq in range(self.mesh.number_of_cell_quadrature_points_in_mesh)
    #     ]
    #     min_eps = min(strain_vals)
    #     max_eps = max(strain_vals)
    #     print("EPS_XY_MIN : {}".format(min_eps))
    #     print("EPS_XY_MAX : {}".format(max_eps))
    #     stress_vals = [
    #         material.mat_data.s1.thermodynamic_forces[nq, 0]
    #         for nq in range(self.mesh.number_of_cell_quadrature_points_in_mesh)
    #     ]
    #     min_sig = min(stress_vals)
    #     max_sig = max(stress_vals)
    #     print("SIG_XX_MIN : {}".format(min_sig))
    #     print("SIG_XX_MAX : {}".format(max_sig))
    #     stress_vals = [
    #         material.mat_data.s1.thermodynamic_forces[nq, 1]
    #         for nq in range(self.mesh.number_of_cell_quadrature_points_in_mesh)
    #     ]
    #     min_sig = min(stress_vals)
    #     max_sig = max(stress_vals)
    #     print("SIG_YY_MIN : {}".format(min_sig))
    #     print("SIG_YY_MAX : {}".format(max_sig))
    #     stress_vals = [
    #         material.mat_data.s1.thermodynamic_forces[nq, 2]
    #         for nq in range(self.mesh.number_of_cell_quadrature_points_in_mesh)
    #     ]
    #     min_sig = min(stress_vals)
    #     max_sig = max(stress_vals)
    #     print("SIG_ZZ_MIN : {}".format(min_sig))
    #     print("SIG_ZZ_MAX : {}".format(max_sig))
    #     stress_vals = [
    #         material.mat_data.s1.thermodynamic_forces[nq, 3]
    #         for nq in range(self.mesh.number_of_cell_quadrature_points_in_mesh)
    #     ]
    #     min_sig = min(stress_vals)
    #     max_sig = max(stress_vals)
    #     print("SIG_XY_MIN : {}".format(min_sig))
    #     print("SIG_XY_MAX : {}".format(max_sig))
