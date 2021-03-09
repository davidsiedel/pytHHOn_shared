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
        # integration_order = 2 * (polynomial_order)
    elif element_type in [ElementType.HHO_EQUAL, ElementType.HDG_EQUAL]:
        integration_order = 2 * (polynomial_order + 1)
        # integration_order = 2 * (polynomial_order)
    elif element_type in [ElementType.HHO_HIGH, ElementType.HDG_HIGH]:
        integration_order = 2 * (polynomial_order + 1)
        # integration_order = 2 * (polynomial_order)
    else:
        raise KeyError("NO")
    print("INTEGRATION ORDER : {}".format(integration_order))
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

    def create_quadrature_points_res_files(self, suffix: str, material: Material):
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
            res_qdp_file.write("STRAIN_TRACE,")
            res_qdp_file.write("HYDRO_STRESS,")
            if material.behaviour_name != "Elasticity":
            # try:
                isv = material.mat_data.s1.internal_state_variables
                for isv_val in range(len(isv[0])):
                    res_qdp_file.write("INTERNAL_STATE_VARIABLE_{},".format(isv_val))
            # except:
            #     pass
            # stored_energies
            # ', '
            # dissipated_energies
            # ', '
            # internal_state_variables
            # '
            # for
            res_qdp_file.write("\n")

    def write_vertex_res_files(self, suffix: str, faces_unknown_vector: ndarray):
        """

        Args:
            suffix:
            faces_unknown_vector:

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
                            # vertex_field_value[u_dir] += self.elements[c].get_cell_field_increment_value(
                            #     point=vertex,
                            #     direction=u_dir,
                            #     field=self.field,
                            #     finite_element=self.finite_element,
                            #     element_unknown_vector=unknown_increment,
                            # )
                            vertex_field_value[u_dir] += self.elements[c].get_cell_field_value(
                                field=self.field,
                                finite_element=self.finite_element,
                                faces_unknown_vector=faces_unknown_vector,
                                point=vertex,
                                direction=u_dir,
                            )
                vertex_field_value = vertex_field_value / self.mesh.vertices_weights_cell[vertex_count]
                for u_dir in range(self.field.field_dimension):
                    res_vtx_file.write("{},".format(vertex_field_value[u_dir]))
                res_vtx_file.write("\n")
        return

    def write_quadrature_points_res_files(self, suffix: str, material: Material, faces_unknown_vector: ndarray):
        """

        Args:
            suffix:
            material:
            faces_unknown_vector:

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
                        # quad_point_field_value = element.get_cell_field_increment_value(
                        #     point=x_q_c,
                        #     direction=u_dir,
                        #     field=self.field,
                        #     finite_element=self.finite_element,
                        #     element_unknown_vector=unknown_increment,
                        # )
                        quad_point_field_value = element.get_cell_field_value(
                            field=self.field,
                            finite_element=self.finite_element,
                            faces_unknown_vector=faces_unknown_vector,
                            point=x_q_c,
                            direction=u_dir,
                        )
                        res_qdp_file.write("{},".format(quad_point_field_value))
                    for g_dir in range(self.field.gradient_dimension):
                        strain_component = material.mat_data.s1.gradients[qp, g_dir]
                        res_qdp_file.write("{},".format(strain_component))
                    for g_dir in range(self.field.gradient_dimension):
                        if self.field.strain_type == StrainType.DISPLACEMENT_TRANSFORMATION_GRADIENT:
                            F = np.zeros((3, 3), dtype=real)
                            F[0, 0] = material.mat_data.s1.gradients[qp, 0]
                            F[1, 1] = material.mat_data.s1.gradients[qp, 1]
                            F[2, 2] = material.mat_data.s1.gradients[qp, 2]
                            F[0, 1] = material.mat_data.s1.gradients[qp, 3]
                            F[1, 0] = material.mat_data.s1.gradients[qp, 4]
                            PK = np.zeros((3, 3),dtype=real)
                            PK[0, 0] = material.mat_data.s1.thermodynamic_forces[qp, 0]
                            PK[1, 1] = material.mat_data.s1.thermodynamic_forces[qp, 1]
                            PK[2, 2] = material.mat_data.s1.thermodynamic_forces[qp, 2]
                            PK[0, 1] = material.mat_data.s1.thermodynamic_forces[qp, 3]
                            PK[1, 0] = material.mat_data.s1.thermodynamic_forces[qp, 4]
                            J = np.linalg.det(F)
                            # F_T_inv = np.linalg.inv(F.T)
                            sig = (1.0 / J) * PK @ F.T
                            sig_vect = np.zeros((5,),dtype=real)
                            sig_vect[0] = sig[0, 0]
                            sig_vect[1] = sig[1, 1]
                            sig_vect[2] = sig[2, 2]
                            sig_vect[3] = sig[0, 1]
                            sig_vect[4] = sig[1, 0]
                            stress_component = sig_vect[g_dir]
                        elif self.field.strain_type == StrainType.DISPLACEMENT_SYMMETRIC_GRADIENT:
                            stress_component = material.mat_data.s1.thermodynamic_forces[qp, g_dir]
                        res_qdp_file.write("{},".format(stress_component))
                    hyrdostatic_pressure = 0.0
                    strain_trace = 0.0
                    if self.field.field_type in [FieldType.DISPLACEMENT_PLANE_STRAIN, FieldType.DISPLACEMENT_PLANE_STRESS]:
                        num_diagonal_components = 3
                    else:
                        num_diagonal_components = self.field.field_dimension
                    for x_dir in range(num_diagonal_components):
                        strain_trace += material.mat_data.s1.gradients[qp, x_dir]
                        hyrdostatic_pressure += material.mat_data.s1.thermodynamic_forces[qp, x_dir]
                    hyrdostatic_pressure = hyrdostatic_pressure/num_diagonal_components
                    res_qdp_file.write("{},".format(strain_trace))
                    res_qdp_file.write("{},".format(hyrdostatic_pressure))
                    # try:
                    #     isv = material.mat_data.s1.internal_state_variables[qp]
                    #     for isv_val in range(len(isv)):
                    #         res_qdp_file.write("{},".format(isv[isv_val]))
                    # except:
                    #     pass
                    if material.behaviour_name != "Elasticity":
                        isv = material.mat_data.s1.internal_state_variables[qp]
                        for isv_val in range(len(isv)):
                            res_qdp_file.write("{},".format(isv[isv_val]))
                    qp += 1
                    res_qdp_file.write("\n")
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

    def solve_newton_2(self, material: Material, verbose: bool, check: bool):
        """

        Args:
            material:

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
        faces_unknown_vector = np.zeros((_constrained_system_size), dtype=real)
        faces_unknown_vector_previous_step = np.zeros((_constrained_system_size), dtype=real)
        residual_values = []
        for time_step_index, time_step in enumerate(self.time_steps):
            # --- SET TEMPERATURE
            material.set_temperature()
            # --- PRINT DATA
            print("----------------------------------------------------------------------------------------------------")
            print("TIME_STEP : {} | LOAD_VALUE : {}".format(time_step_index, time_step))
            # --- WRITE RES FILES
            file_suffix = "{}".format(time_step_index).zfill(6)
            self.create_vertex_res_files(file_suffix)
            self.create_quadrature_points_res_files(file_suffix, material)
            # correction = np.zeros((_constrained_system_size),dtype=real)
            # --- RESET DISPLACEMENT
            reset_displacement = False
            if reset_displacement:
                faces_unknown_vector = np.zeros((_constrained_system_size),dtype=real)
                for element in self.elements:
                    element.cell_unknown_vector = np.zeros((_dx * _cl, ),dtype=real)
            for iteration in range(self.number_of_iterations):
                if check:
                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    noise_tolerance = 1.0e-6
                    noise = noise_tolerance * np.max(np.abs(np.copy(faces_unknown_vector[:_system_size])))
                    numerical_forces_bs_0 = [
                        np.zeros((element.element_size, element.element_size), dtype=real) for element in self.elements
                    ]
                    numerical_forces_as_0 = [
                        np.zeros((element.element_size, element.element_size), dtype=real) for element in self.elements
                    ]
                    numerical_unknowns_0 = [
                        np.zeros((element.element_size, element.element_size), dtype=real) for element in self.elements
                    ]
                    for epsi_dir in range(self.elements[0].element_size):
                        _qp = 0
                        for _element_index, element in enumerate(self.elements):
                            # --- GET ELEMENT UNKNOWN INCREMENT
                            element_unknown_vector = element.get_element_unknown_vector(
                                self.field, self.finite_element, faces_unknown_vector
                            )
                            numerical_unknowns_0[_element_index][:, epsi_dir] += np.copy(element_unknown_vector)
                            numerical_unknowns_0[_element_index][epsi_dir, epsi_dir] += noise
                            for _qc in range(len(element.cell.quadrature_weights)):
                                _w_q_c = element.cell.quadrature_weights[_qc]
                                _x_q_c = element.cell.quadrature_points[:, _qc]
                                transformation_gradient_0 = element.gradients_operators[_qc] @ numerical_unknowns_0[_element_index][:, epsi_dir]
                                material.mat_data.s1.gradients[_qp] = transformation_gradient_0
                                integ_res = mgis_bv.integrate(
                                    material.mat_data, material.integration_type, 0.0, _qp, (_qp + 1)
                                )
                                numerical_forces_bs_0[_element_index][:, epsi_dir] += _w_q_c * (
                                        element.gradients_operators[_qc].T @ material.mat_data.s1.thermodynamic_forces[_qp]
                                )
                                _qp += 1
                        mgis_bv.revert(material.mat_data)
                    for _element_index, element in enumerate(self.elements):
                        numerical_forces_as_0[_element_index] += numerical_forces_bs_0[_element_index]
                        for epsi_dir in range(self.elements[0].element_size):
                            numerical_forces_as_0[_element_index][:, epsi_dir] += (
                                    material.stabilization_parameter
                                    * element.stabilization_operator
                                    @ numerical_unknowns_0[_element_index][:, epsi_dir]
                            )
                    numerical_forces_bs_1 = [
                        np.zeros((element.element_size, element.element_size), dtype=real) for element in self.elements
                    ]
                    numerical_forces_as_1 = [
                        np.zeros((element.element_size, element.element_size), dtype=real) for element in self.elements
                    ]
                    numerical_unknowns_1 = [
                        np.zeros((element.element_size, element.element_size), dtype=real) for element in self.elements
                    ]
                    for epsi_dir in range(self.elements[0].element_size):
                        _qp = 0
                        for _element_index, element in enumerate(self.elements):
                            # --- GET ELEMENT UNKNOWN INCREMENT
                            element_unknown_increment = element.get_element_unknown_vector(
                                self.field, self.finite_element, faces_unknown_vector
                            )
                            numerical_unknowns_1[_element_index][:, epsi_dir] += np.copy(element_unknown_increment)
                            numerical_unknowns_1[_element_index][epsi_dir, epsi_dir] -= noise
                            for _qc in range(len(element.cell.quadrature_weights)):
                                _w_q_c = element.cell.quadrature_weights[_qc]
                                _x_q_c = element.cell.quadrature_points[:, _qc]
                                transformation_gradient_1 = element.gradients_operators[_qc] @ numerical_unknowns_1[_element_index][:, epsi_dir]
                                material.mat_data.s1.gradients[_qp] = transformation_gradient_1
                                integ_res = mgis_bv.integrate(
                                    material.mat_data, material.integration_type, 0.0, _qp, (_qp + 1)
                                )
                                numerical_forces_bs_1[_element_index][:, epsi_dir] += _w_q_c * (
                                        element.gradients_operators[_qc].T @ material.mat_data.s1.thermodynamic_forces[_qp]
                                )
                                _qp += 1
                        mgis_bv.revert(material.mat_data)
                    for _element_index, element in enumerate(self.elements):
                        numerical_forces_as_1[_element_index] = +numerical_forces_bs_1[_element_index]
                        for epsi_dir in range(self.elements[0].element_size):
                            numerical_forces_as_1[_element_index][:, epsi_dir] += (
                                    material.stabilization_parameter
                                    * element.stabilization_operator
                                    @ numerical_unknowns_1[_element_index][:, epsi_dir]
                            )
                    # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # --------------------------------------------------------------------------------------------------
                # SET SYSTEM MATRIX AND VECTOR
                # --------------------------------------------------------------------------------------------------
                tangent_matrix = np.zeros((_constrained_system_size, _constrained_system_size), dtype=real)
                residual = np.zeros((_constrained_system_size), dtype=real)
                # --------------------------------------------------------------------------------------------------
                # SET TIME INCREMENT
                # --------------------------------------------------------------------------------------------------
                if time_step_index == 0:
                    _dt = time_step
                else:
                    _dt = time_step - self.time_steps[time_step_index - 1]
                    _dt = np.float64(_dt)
                # _dt = 0.0
                # _dt = np.float64(0.0)
                # --------------------------------------------------------------------------------------------------
                # FOR ELEMENT LOOP
                # --------------------------------------------------------------------------------------------------
                _qp = 0
                stab_coef = 1.0
                stab_inf = np.inf
                iter_face_constraint = 0
                if False:
                    for element in self.elements:
                        for _qc in range(len(element.cell.quadrature_weights)):
                            _w_q_c = element.cell.quadrature_weights[_qc]
                            _x_q_c = element.cell.quadrature_points[:, _qc]
                            # --- COMPUTE STRAINS AND SET THEM IN THE BEHAVIOUR LAW
                            transformation_gradient = element.get_transformation_gradient(
                                self.field, self.finite_element, faces_unknown_vector, _qc
                            )
                            material.mat_data.s1.gradients[_qp] = transformation_gradient
                            # --- INTEGRATE BEHAVIOUR LAW
                            integ_res = mgis_bv.integrate(material.mat_data, material.integration_type, _dt, _qp, (_qp + 1))
                            eigen_values = np.linalg.eig(material.mat_data.K[_qp])[0]
                            stab_elem = np.min(eigen_values)/material.stabilization_parameter
                            if stab_elem < stab_inf and stab_elem > 1.0:
                                stab_inf = stab_elem
                                stab_coef = stab_elem
                            _qp += 1
                    mgis_bv.revert(material.mat_data)
                _qp = 0
                for _element_index, element in enumerate(self.elements):
                    _nf = len(element.faces)
                    _c0_c = _dx * _cl
                    # --- INITIALIZE MATRIX AND VECTORS
                    element_stiffness_matrix = np.zeros((element.element_size, element.element_size), dtype=real)
                    element_internal_forces = np.zeros((element.element_size,), dtype=real)
                    element_external_forces = np.zeros((element.element_size,), dtype=real)
                    # --- RUN OVER EACH QUADRATURE POINT
                    for _qc in range(len(element.cell.quadrature_weights)):
                        _w_q_c = element.cell.quadrature_weights[_qc]
                        _x_q_c = element.cell.quadrature_points[:, _qc]
                        # --- COMPUTE STRAINS AND SET THEM IN THE BEHAVIOUR LAW
                        transformation_gradient = element.get_transformation_gradient(
                            self.field, self.finite_element, faces_unknown_vector, _qc
                        )
                        material.mat_data.s1.gradients[_qp] = transformation_gradient
                        # --- INTEGRATE BEHAVIOUR LAW
                        integ_res = mgis_bv.integrate(material.mat_data, material.integration_type, _dt, _qp, (_qp + 1))
                        # stored_energies', 'dissipated_energies', 'internal_state_variables
                        # print(material.mat_data.s1.internal_state_variables)
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
                    if verbose:
                        print("ELEM : {} | INTERNAL_FORCES_BEFORE STAB : \n {}".format(_element_index, element_internal_forces))
                    # --- STAB PARAMETER CHANGE
                    stab_param = stab_coef * material.stabilization_parameter
                    if check and iteration > 0:
                        numerical_element_stiffness_matrix = (
                            numerical_forces_bs_0[_element_index] - numerical_forces_bs_1[_element_index]
                        ) / (2.0 * noise)
                        check_num = np.max(
                            np.abs(numerical_element_stiffness_matrix - element_stiffness_matrix)
                        ) / np.max(np.abs(element_stiffness_matrix))
                        if check_num > noise_tolerance:
                            print(
                                "AT ELEM : {} | CHECK_BEFORE_STAB : {}".format(str(_element_index).zfill(6), check_num)
                            )
                    # --- ADDING STABILIZATION CONTRIBUTION AT THE ELEMENT LEVEL
                    element_stiffness_matrix += stab_param * element.stabilization_operator
                    # element_stiffness_matrix -= stab_param * element.stabilization_operator
                    # --- ADDING STABILIZATION CONTRIBUTION AT THE ELEMENT LEVEL
                    element_internal_forces += (
                        stab_param
                        * element.stabilization_operator
                        @ element.get_element_unknown_vector(self.field, self.finite_element, faces_unknown_vector)
                    )
                    # element_internal_forces -= (
                    #     stab_param
                    #     * element.stabilization_operator
                    #     @ element.get_element_unknown_vector(self.field, self.finite_element, faces_unknown_vector)
                    # )
                    if check and iteration > 0:
                        numerical_element_stiffness_matrix = (
                            numerical_forces_as_0[_element_index] - numerical_forces_as_1[_element_index]
                        ) / (2.0 * noise)
                        check_num = np.max(
                            np.abs(numerical_element_stiffness_matrix - element_stiffness_matrix)
                        ) / np.max(np.abs(element_stiffness_matrix))
                        if check_num > noise_tolerance:
                            print(
                                "AT ELEM : {} | CHECK_BEFORE_STAB : {}".format(str(_element_index).zfill(6), check_num)
                            )
                    if verbose:
                        print("ELEM : {} | INTERNAL_FORCES_AFTER STAB : \n {}".format(_element_index, element_internal_forces))
                        _iv0 = _element_index * len(element.cell.quadrature_weights)
                        _iv1 = (_element_index + 1) * len(element.cell.quadrature_weights)
                        print("ELEM : {} | DISPLACEMENT S0 : \n {}".format(_element_index, element.get_element_unknown_vector(self.field, self.finite_element, faces_unknown_vector_previous_step)))
                        print("ELEM : {} | GRADIENTS S0 : \n {}".format(_element_index, material.mat_data.s0.gradients[_iv0:_iv1]))
                        print("ELEM : {} | DISPLACEMENT S1 : \n {}".format(_element_index, element.get_element_unknown_vector(self.field, self.finite_element, faces_unknown_vector)))
                        print("ELEM : {} | GRADIENTS S1 : \n {}".format(_element_index, material.mat_data.s1.gradients[_iv0:_iv1]))
                        if not material.behaviour_name == "Elasticity":
                            print("ELEM : {} | INTERNAL_STATE_VARIABLES S0 : \n {}".format(_element_index, material.mat_data.s0.internal_state_variables[_iv0:_iv1]))
                            print("ELEM : {} | INTERNAL_STATE_VARIABLES S1 : \n {}".format(_element_index, material.mat_data.s1.internal_state_variables[_iv0:_iv1]))
                    # --- BOUNDARY CONDITIONS
                    for boundary_condition in self.boundary_conditions:
                        # --- DISPLACEMENT CONDITIONS
                        if boundary_condition.boundary_type == BoundaryType.DISPLACEMENT:
                            for f_local, f_global in enumerate(element.faces_indices):
                                if (
                                    f_global
                                    in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]
                                ):
                                    _l0 = _system_size + iter_face_constraint * _fk
                                    _l1 = _system_size + (iter_face_constraint + 1) * _fk
                                    _c0 = _cl * _dx + (f_local * _dx * _fk) + boundary_condition.direction * _fk
                                    _c1 = _cl * _dx + (f_local * _dx * _fk) + (boundary_condition.direction + 1) * _fk
                                    _r0 = f_global * _fk * _dx + _fk * boundary_condition.direction
                                    _r1 = f_global * _fk * _dx + _fk * (boundary_condition.direction + 1)
                                    # -------
                                    # face_displacement = element_unknown_increment[_r0:_r1]
                                    # face_displacement = np.copy(element_unknown_increment[_c0:_c1])
                                    # face_displacement = element_unknown_increment[_c0:_c1]
                                    face_lagrange = faces_unknown_vector[_l0:_l1]
                                    face_displacement = faces_unknown_vector[_r0:_r1]
                                    _m_psi_psi_face = np.zeros((_fk, _fk), dtype=real)
                                    _v_face_imposed_displacement = np.zeros((_fk,), dtype=real)
                                    for _qf in range(len(element.faces[f_local].quadrature_weights)):
                                        _h_f = element.faces[f_local].shape.diameter
                                        _x_f = element.faces[f_local].shape.centroid
                                        _x_q_f = element.faces[f_local].quadrature_points[:, _qf]
                                        _w_q_f = element.faces[f_local].quadrature_weights[_qf]
                                        _s_q_f = (element.faces[f_local].mapping_matrix @ _x_q_f)[:-1]
                                        _s_f = (element.faces[f_local].mapping_matrix @ _x_f)[:-1]
                                        # v = self.finite_element.face_basis_k.evaluate_function(
                                        #     _x_q_f,
                                        #     element.faces[f_local].shape.centroid,
                                        #     element.faces[f_local].shape.diameter,
                                        # )
                                        v = self.finite_element.face_basis_k.evaluate_function(
                                            _s_q_f,
                                            _s_f,
                                            _h_f,
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
                                    if debug_mode == 0:
                                        print("FACE MASS MATRIX IN DIRICHLET BOUND COND :")
                                        print("{}".format(np.linalg.cond(_m_psi_psi_face)))
                                    imposed_face_displacement = _m_psi_psi_face_inv @ _v_face_imposed_displacement
                                    face_displacement_difference = face_displacement - imposed_face_displacement
                                    # --- LAGRANGE INTERNAL FORCES PART
                                    element_internal_forces[_c0:_c1] += (
                                        normalization_lagrange_coefficient * face_lagrange
                                    )
                                    # --- LAGRANGE MULTIPLIERS PART
                                    # residual[_l0:_l1] += (
                                    #     normalization_lagrange_coefficient * face_displacement_difference
                                    # )
                                    residual[_l0:_l1] -= (
                                        normalization_lagrange_coefficient * face_displacement_difference
                                    )
                                    # --- LAGRANGE MATRIX PART
                                    tangent_matrix[_l0:_l1, _r0:_r1] += normalization_lagrange_coefficient * np.eye(_fk, dtype=real)
                                    tangent_matrix[_r0:_r1, _l0:_l1] += normalization_lagrange_coefficient * np.eye(_fk, dtype=real)
                                    # --- SET EXTERNAL FORCES COEFFICIENT
                                    lagrange_external_forces = (
                                        normalization_lagrange_coefficient * imposed_face_displacement
                                    )
                                    if np.max(np.abs(lagrange_external_forces)) > external_forces_coefficient:
                                        external_forces_coefficient = np.max(np.abs(lagrange_external_forces))
                                    iter_face_constraint += 1
                        elif boundary_condition.boundary_type == BoundaryType.PRESSURE:
                            for f_local, f_global in enumerate(element.faces_indices):
                                if (
                                    f_global
                                    in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]
                                ):
                                    for qf in range(len(element.faces[f_local].quadrature_weights)):
                                        _h_f = element.faces[f_local].shape.diameter
                                        _x_f = element.faces[f_local].shape.centroid
                                        _x_q_f = element.faces[f_local].quadrature_points[:, _qf]
                                        _w_q_f = element.faces[f_local].quadrature_weights[_qf]
                                        _s_q_f = (element.faces[f_local].mapping_matrix @ _x_q_f)[:-1]
                                        _s_f = (element.faces[f_local].mapping_matrix @ _x_f)[:-1]
                                        # v = self.finite_element.face_basis_k.evaluate_function(
                                        #     _x_q_f,
                                        #     element.faces[f_local].shape.centroid,
                                        #     element.faces[f_local].shape.diameter,
                                        # )
                                        v = self.finite_element.face_basis_k.evaluate_function(
                                            _s_q_f,
                                            _s_f,
                                            _h_f,
                                        )
                                        # _x_q_f = element.faces[f_local].quadrature_points[:, qf]
                                        # _w_q_f = element.faces[f_local].quadrature_weights[qf]
                                        # v = self.finite_element.face_basis_k.evaluate_function(
                                        #     _x_q_f,
                                        #     element.faces[f_local].shape.centroid,
                                        #     element.faces[f_local].shape.diameter,
                                        # )
                                        vf = _w_q_f * v * boundary_condition.function(time_step, _x_q_f)
                                        _c0 = _dx * _cl + f_local * _dx * _fk + boundary_condition.direction * _fk
                                        _c1 = _dx * _cl + f_local * _dx * _fk + (boundary_condition.direction + 1) * _fk
                                        # _r0 = f_global * _fk * _dx + _fk * boundary_condition.direction
                                        # _r1 = f_global * _fk * _dx + _fk * (boundary_condition.direction + 1)
                                        element_external_forces[_c0:_c1] += vf
                    # --- COMPUTE RESIDUAL AFTER VOLUMETRIC CONTRIBUTION
                    element_residual = element_internal_forces - element_external_forces
                    if verbose:
                        print("ELEM : {} | INTERNAL_FORCES_END : \n {}".format(_element_index,element_internal_forces))
                        print("ELEM : {} | ELEMENT_RESIDUAL : \n {}".format(_element_index, element_residual))
                    # --------------------------------------------------------------------------------------------------
                    # CONDENSATION
                    # --------------------------------------------------------------------------------------------------
                    m_cell_cell = element_stiffness_matrix[:_c0_c, :_c0_c]
                    m_cell_faces = element_stiffness_matrix[:_c0_c, _c0_c:]
                    m_faces_cell = element_stiffness_matrix[_c0_c:, :_c0_c]
                    m_faces_faces = element_stiffness_matrix[_c0_c:, _c0_c:]
                    # v_cell = element_residual[:_c0_c]
                    # v_faces = element_residual[_c0_c:]
                    v_cell = -element_residual[:_c0_c]
                    v_faces = -element_residual[_c0_c:]
                    _condtest = np.linalg.cond(m_cell_cell)
                    m_cell_cell_inv = np.linalg.inv(m_cell_cell)
                    if debug_mode == 0:
                        print("TANGENT MATRIX IN CONDENSATION COND :")
                        print("{}".format(np.linalg.cond(m_cell_cell)))
                    # ge = m_faces_cell @ m_cell_cell_inv
                    # gd = (m_faces_cell @ m_cell_cell_inv) @ m_cell_faces
                    K_cond = m_faces_faces - ((m_faces_cell @ m_cell_cell_inv) @ m_cell_faces)
                    R_cond = v_faces - (m_faces_cell @ m_cell_cell_inv) @ v_cell
                    # --- SET CONDENSATION/DECONDENSATION MATRICES
                    element.m_cell_cell_inv = m_cell_cell_inv
                    element.m_cell_faces = m_cell_faces
                    element.v_cell = v_cell
                    # --- ASSEMBLY
                    for _i_local, _i_global in enumerate(element.faces_indices):
                        _rg0 = _i_global * (_fk * _dx)
                        _rg1 = (_i_global + 1) * (_fk * _dx)
                        _re0 = _i_local * (_fk * _dx)
                        _re1 = (_i_local + 1) * (_fk * _dx)
                        residual[_rg0:_rg1] += R_cond[_re0:_re1]
                        for _j_local, _j_global in enumerate(element.faces_indices):
                            _cg0 = _j_global * (_fk * _dx)
                            _cg1 = (_j_global + 1) * (_fk * _dx)
                            _ce0 = _j_local * (_fk * _dx)
                            _ce1 = (_j_local + 1) * (_fk * _dx)
                            tangent_matrix[_rg0:_rg1, _cg0:_cg1] += K_cond[_re0:_re1, _ce0:_ce1]
                    # --- SET EXTERNAL FORCES COEFFICIENT
                    if np.max(np.abs(element_external_forces)) > external_forces_coefficient:
                        external_forces_coefficient = np.max(np.abs(element_external_forces))
                # --------------------------------------------------------------------------------------------------
                # RESIDUAL EVALUATION
                # --------------------------------------------------------------------------------------------------
                if external_forces_coefficient == 0.0:
                    external_forces_coefficient = 1.0
                residual_evaluation = np.max(np.abs(residual)) / external_forces_coefficient
                print("ITER : {} | RES_MAX : {}".format(str(iteration).zfill(4), residual_evaluation))
                # print("ITER : {} | =====================================================================".format(str(iteration).zfill(4)))
                # residual_values.append(residual_evaluation)
                # if residual_evaluation < self.tolerance:
                if residual_evaluation < self.tolerance or iteration == self.number_of_iterations - 1:
                    # ----------------------------------------------------------------------------------------------
                    # UPDATE INTERNAL VARIABLES
                    # ----------------------------------------------------------------------------------------------
                    mgis_bv.update(material.mat_data)
                    print("ITERATIONS : {}".format(iteration + 1))
                    self.write_vertex_res_files(file_suffix, faces_unknown_vector)
                    self.write_quadrature_points_res_files(file_suffix, material, faces_unknown_vector)
                    faces_unknown_vector_previous_step += faces_unknown_vector
                    residual_values.append(residual_evaluation)
                    break
                else:
                    # ----------------------------------------------------------------------------------------------
                    # SOLVE SYSTEM
                    # ----------------------------------------------------------------------------------------------
                    if debug_mode == 0:
                        print("K GLOBAL COND")
                        print("{}".format(np.linalg.cond(tangent_matrix)))
                    # sparse_global_matrix = csr_matrix(-tangent_matrix)
                    sparse_global_matrix = csr_matrix(tangent_matrix)
                    correction = spsolve(sparse_global_matrix, residual)
                    # print("CORRECTION :")
                    # print(correction)
                    # correction = np.linalg.solve(-tangent_matrix, residual)
                    faces_unknown_vector += correction
                    if verbose:
                        print("R_K_RES : \n {}".format(residual - tangent_matrix @ correction))
                        print("R_K_RES_NORMALIZED : \n {}".format((residual - tangent_matrix @ correction)/external_forces_coefficient))
                    # ----------------------------------------------------------------------------------------------
                    # DECONDENSATION
                    # ----------------------------------------------------------------------------------------------
                    for element in self.elements:
                        _nf = len(element.faces)
                        # _c0_c = _dx * _cl
                        # --- DECONDENSATION - GET ELEMENT UNKNOWN CORRECTION
                        face_correction = np.zeros((_nf * _fk * _dx), dtype=real)
                        # element_correction = np.zeros((element.element_size,))
                        for _i_local, _i_global in enumerate(element.faces_indices):
                            _c0_fg = _i_global * (_fk * _dx)
                            _c1_fg = (_i_global + 1) * (_fk * _dx)
                            # _c0_fl = (_dx * _cl) + _i_local * (_fk * _dx)
                            # _c1_fl = (_dx * _cl) + (_i_local + 1) * (_fk * _dx)
                            _c0_fl = _i_local * (_fk * _dx)
                            _c1_fl = (_i_local + 1) * (_fk * _dx)
                            # element_correction[_c0_fl:_c1_fl] += correction[_c0_fg:_c1_fg]
                            face_correction[_c0_fl:_c1_fl] += correction[_c0_fg:_c1_fg]
                        # element_correction[:_c0_c] = element.m_cell_cell_inv @ (
                        #     element.v_cell - element.m_cell_faces @ element_correction[_c0_c:]
                        # )
                        cell_correction = element.m_cell_cell_inv @ (
                            element.v_cell - element.m_cell_faces @ face_correction
                        )
                        # --- ADDING CORRECTION TO CURRENT DISPLACEMENT
                        element.cell_unknown_vector += cell_correction
        plt.plot(range(len(residual_values)), residual_values, label="residual")
        plt.plot(range(len(residual_values)), [self.tolerance for local_i in range(len(residual_values))], label="tolerance")
        plt.ylabel("resiudal after convergence")
        plt.xlabel("time step index")
        plt.title("evolution of the residual")
        plt.legend()
        plt.gcf().subplots_adjust(left=0.15)
        plt.show()
        return

    def solve_newton_exact(self, material: Material, verbose: bool, check: bool):
        """

        Args:
            material:

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
        faces_unknown_vector = np.zeros((_constrained_system_size), dtype=real)
        faces_unknown_vector_previous_step = np.zeros((_constrained_system_size), dtype=real)
        residual_values = []
        correction = np.zeros((_constrained_system_size), dtype=real)
        for time_step_index, time_step in enumerate(self.time_steps):
            # --- SET TEMPERATURE
            material.set_temperature()
            # --- PRINT DATA
            print("----------------------------------------------------------------------------------------------------")
            print("TIME_STEP : {} | LOAD_VALUE : {}".format(time_step_index, time_step))
            # --- WRITE RES FILES
            file_suffix = "{}".format(time_step_index).zfill(6)
            self.create_vertex_res_files(file_suffix)
            self.create_quadrature_points_res_files(file_suffix, material)
            # correction = np.zeros((_constrained_system_size),dtype=real)
            # --- RESET DISPLACEMENT
            reset_displacement = False
            if reset_displacement:
                faces_unknown_vector = np.zeros((_constrained_system_size),dtype=real)
                for element in self.elements:
                    element.cell_unknown_vector = np.zeros((_dx * _cl, ),dtype=real)
            for iteration in range(self.number_of_iterations):
                if check:
                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    noise_tolerance = 1.0e-6
                    noise = noise_tolerance * np.max(np.abs(np.copy(faces_unknown_vector[:_system_size])))
                    numerical_forces_bs_0 = [
                        np.zeros((element.element_size, element.element_size), dtype=real) for element in self.elements
                    ]
                    numerical_forces_as_0 = [
                        np.zeros((element.element_size, element.element_size), dtype=real) for element in self.elements
                    ]
                    numerical_unknowns_0 = [
                        np.zeros((element.element_size, element.element_size), dtype=real) for element in self.elements
                    ]
                    for epsi_dir in range(self.elements[0].element_size):
                        _qp = 0
                        for _element_index, element in enumerate(self.elements):
                            # --- GET ELEMENT UNKNOWN INCREMENT
                            element_unknown_vector = element.get_element_unknown_vector(
                                self.field, self.finite_element, faces_unknown_vector
                            )
                            numerical_unknowns_0[_element_index][:, epsi_dir] += np.copy(element_unknown_vector)
                            numerical_unknowns_0[_element_index][epsi_dir, epsi_dir] += noise
                            for _qc in range(len(element.cell.quadrature_weights)):
                                _w_q_c = element.cell.quadrature_weights[_qc]
                                _x_q_c = element.cell.quadrature_points[:, _qc]
                                transformation_gradient_0 = element.gradients_operators[_qc] @ numerical_unknowns_0[_element_index][:, epsi_dir]
                                material.mat_data.s1.gradients[_qp] = transformation_gradient_0
                                integ_res = mgis_bv.integrate(
                                    material.mat_data, material.integration_type, 0.0, _qp, (_qp + 1)
                                )
                                numerical_forces_bs_0[_element_index][:, epsi_dir] += _w_q_c * (
                                        element.gradients_operators[_qc].T @ material.mat_data.s1.thermodynamic_forces[_qp]
                                )
                                _qp += 1
                        mgis_bv.revert(material.mat_data)
                    for _element_index, element in enumerate(self.elements):
                        numerical_forces_as_0[_element_index] += numerical_forces_bs_0[_element_index]
                        for epsi_dir in range(self.elements[0].element_size):
                            numerical_forces_as_0[_element_index][:, epsi_dir] += (
                                    material.stabilization_parameter
                                    * element.stabilization_operator
                                    @ numerical_unknowns_0[_element_index][:, epsi_dir]
                            )
                    numerical_forces_bs_1 = [
                        np.zeros((element.element_size, element.element_size), dtype=real) for element in self.elements
                    ]
                    numerical_forces_as_1 = [
                        np.zeros((element.element_size, element.element_size), dtype=real) for element in self.elements
                    ]
                    numerical_unknowns_1 = [
                        np.zeros((element.element_size, element.element_size), dtype=real) for element in self.elements
                    ]
                    for epsi_dir in range(self.elements[0].element_size):
                        _qp = 0
                        for _element_index, element in enumerate(self.elements):
                            # --- GET ELEMENT UNKNOWN INCREMENT
                            element_unknown_increment = element.get_element_unknown_vector(
                                self.field, self.finite_element, faces_unknown_vector
                            )
                            numerical_unknowns_1[_element_index][:, epsi_dir] += np.copy(element_unknown_increment)
                            numerical_unknowns_1[_element_index][epsi_dir, epsi_dir] -= noise
                            for _qc in range(len(element.cell.quadrature_weights)):
                                _w_q_c = element.cell.quadrature_weights[_qc]
                                _x_q_c = element.cell.quadrature_points[:, _qc]
                                transformation_gradient_1 = element.gradients_operators[_qc] @ numerical_unknowns_1[_element_index][:, epsi_dir]
                                material.mat_data.s1.gradients[_qp] = transformation_gradient_1
                                integ_res = mgis_bv.integrate(
                                    material.mat_data, material.integration_type, 0.0, _qp, (_qp + 1)
                                )
                                numerical_forces_bs_1[_element_index][:, epsi_dir] += _w_q_c * (
                                        element.gradients_operators[_qc].T @ material.mat_data.s1.thermodynamic_forces[_qp]
                                )
                                _qp += 1
                        mgis_bv.revert(material.mat_data)
                    for _element_index, element in enumerate(self.elements):
                        numerical_forces_as_1[_element_index] = +numerical_forces_bs_1[_element_index]
                        for epsi_dir in range(self.elements[0].element_size):
                            numerical_forces_as_1[_element_index][:, epsi_dir] += (
                                    material.stabilization_parameter
                                    * element.stabilization_operator
                                    @ numerical_unknowns_1[_element_index][:, epsi_dir]
                            )
                    # CHECK ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # --------------------------------------------------------------------------------------------------
                # SET SYSTEM MATRIX AND VECTOR
                # --------------------------------------------------------------------------------------------------
                tangent_matrix = np.zeros((_constrained_system_size, _constrained_system_size), dtype=real)
                residual = np.zeros((_constrained_system_size), dtype=real)
                # --------------------------------------------------------------------------------------------------
                # SET TIME INCREMENT
                # --------------------------------------------------------------------------------------------------
                if time_step_index == 0:
                    _dt = time_step
                else:
                    _dt = time_step - self.time_steps[time_step_index - 1]
                    _dt = np.float64(_dt)
                # --------------------------------------------------------------------------------------------------
                # FOR ELEMENT LOOP
                # --------------------------------------------------------------------------------------------------
                _qp = 0
                stab_coef = 1.0
                stab_inf = np.inf
                iter_face_constraint = 0
                if False:
                    for element in self.elements:
                        for _qc in range(len(element.cell.quadrature_weights)):
                            _w_q_c = element.cell.quadrature_weights[_qc]
                            _x_q_c = element.cell.quadrature_points[:, _qc]
                            # --- COMPUTE STRAINS AND SET THEM IN THE BEHAVIOUR LAW
                            transformation_gradient = element.get_transformation_gradient(
                                self.field, self.finite_element, faces_unknown_vector, _qc
                            )
                            material.mat_data.s1.gradients[_qp] = transformation_gradient
                            # --- INTEGRATE BEHAVIOUR LAW
                            integ_res = mgis_bv.integrate(material.mat_data, material.integration_type, _dt, _qp, (_qp + 1))
                            eigen_values = np.linalg.eig(material.mat_data.K[_qp])[0]
                            stab_elem = np.min(eigen_values)/material.stabilization_parameter
                            if stab_elem < stab_inf and stab_elem > 1.0:
                                stab_inf = stab_elem
                                stab_coef = stab_elem
                            _qp += 1
                        mgis_bv.revert(material.mat_data)
                _qp = 0
                for _element_index, element in enumerate(self.elements):
                    print("ELEMENT : {}".format(_element_index))
                    # --- GET FACES LOCAL CORRECTION
                    local_faces_correction = np.zeros((len(element.faces) * _dx * _fk), dtype=real)
                    for _f_local, _f_global in enumerate(element.faces_indices):
                        local_faces_correction[_f_local * _dx * _fk:(_f_local + 1) * _dx * _fk] += (
                            correction[_f_global * _dx * _fk:(_f_global + 1) * _dx * _fk]
                        )
                    # --- LOCAL NEWTON
                    number_of_local_iterations = 20
                    # local_tolerance = self.tolerance
                    local_tolerance = 1.e-4
                    R_cell_value_previous = np.inf
                    for _local_iteration in range(number_of_local_iterations):
                        _nf = len(element.faces)
                        _c0_c = _dx * _cl
                        # --- INITIALIZE MATRIX AND VECTORS
                        element_stiffness_matrix = np.zeros((element.element_size, element.element_size), dtype=real)
                        element_internal_forces = np.zeros((element.element_size,), dtype=real)
                        element_external_forces = np.zeros((element.element_size,), dtype=real)
                        # --- RUN OVER EACH QUADRATURE POINT
                        _local_qp = 0
                        for _qc in range(len(element.cell.quadrature_weights)):
                            _w_q_c = element.cell.quadrature_weights[_qc]
                            _x_q_c = element.cell.quadrature_points[:, _qc]
                            # --- COMPUTE STRAINS AND SET THEM IN THE BEHAVIOUR LAW
                            transformation_gradient = element.get_transformation_gradient(
                                self.field, self.finite_element, faces_unknown_vector, _qc
                            )
                            material.mat_data.s1.gradients[_qp] = transformation_gradient
                            # --- INTEGRATE BEHAVIOUR LAW
                            integ_res = mgis_bv.integrate(material.mat_data, material.integration_type, _dt, _qp, (_qp + 1))
                            # stored_energies', 'dissipated_energies', 'internal_state_variables
                            # print(material.mat_data.s1.internal_state_variables)
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
                            _local_qp += 1
                        if verbose:
                            print("ELEM : {} | INTERNAL_FORCES_BEFORE STAB : \n {}".format(_element_index, element_internal_forces))
                        # --- STAB PARAMETER CHANGE
                        stab_param = stab_coef * material.stabilization_parameter
                        if check and iteration > 0:
                            numerical_element_stiffness_matrix = (
                                numerical_forces_bs_0[_element_index] - numerical_forces_bs_1[_element_index]
                            ) / (2.0 * noise)
                            check_num = np.max(
                                np.abs(numerical_element_stiffness_matrix - element_stiffness_matrix)
                            ) / np.max(np.abs(element_stiffness_matrix))
                            if check_num > noise_tolerance:
                                print(
                                    "AT ELEM : {} | CHECK_BEFORE_STAB : {}".format(str(_element_index).zfill(6), check_num)
                                )
                        # --- ADDING STABILIZATION CONTRIBUTION AT THE ELEMENT LEVEL
                        element_stiffness_matrix += stab_param * element.stabilization_operator
                        # element_stiffness_matrix -= stab_param * element.stabilization_operator
                        # --- ADDING STABILIZATION CONTRIBUTION AT THE ELEMENT LEVEL
                        element_internal_forces += (
                            stab_param
                            * element.stabilization_operator
                            @ element.get_element_unknown_vector(self.field, self.finite_element, faces_unknown_vector)
                        )
                        # --- MATRIX VECTOR DECOMPOSITION
                        K_cc = element_stiffness_matrix[:_c0_c,:_c0_c]
                        K_cf = element_stiffness_matrix[:_c0_c,_c0_c:]
                        # local_external_forces_coefficient = np.max(np.abs(element_external_forces))
                        element_residual = element_internal_forces - element_external_forces
                        # if local_external_forces_coefficient == 0.:
                        #     local_external_forces_coefficient = 1.0
                            # local_external_forces_coefficient = normalization_lagrange_coefficient
                            # local_external_forces_coefficient = 1./element.cell.shape.volume
                        R_cc = element_residual[:_c0_c]
                        if iteration == 0:
                            local_external_forces_coefficient = 1.0
                        if _local_iteration == 0 and iteration > 0:
                            local_external_forces_coefficient = np.max(np.abs(R_cc))
                            # local_external_forces_coefficient = 1.
                            # if local_external_forces_coefficient == 0.0 or local_external_forces_coefficient < local_tolerance:
                            if local_external_forces_coefficient == 0.0:
                                local_external_forces_coefficient = 1.0
                        # R_cell = R_cc - K_cf @ local_faces_correction
                        # print("LOCAL_EXTERNAL_FORCES_COEF : {}".format(local_external_forces_coefficient))
                        R_cell = R_cc
                        local_residual_evaluation = np.max(np.abs(R_cell) / local_external_forces_coefficient)
                        if local_residual_evaluation > R_cell_value_previous:
                            print("!!!! RESIUDAL INCREASE :\n {}".format(element.cell.quadrature_points))
                        R_cell_value_previous = local_residual_evaluation
                        print("LOCAL ITER : {} | RES MAX : {}".format(_local_iteration, local_residual_evaluation))
                        if local_residual_evaluation < local_tolerance:
                            break
                        else:
                            cell_correction = np.linalg.solve(-K_cc, R_cell)
                            element.cell_unknown_vector += cell_correction
                            if _local_iteration != number_of_local_iterations - 1:
                                _qp -= _local_qp
                    # element_internal_forces -= (
                    #     stab_param
                    #     * element.stabilization_operator
                    #     @ element.get_element_unknown_vector(self.field, self.finite_element, faces_unknown_vector)
                    # )
                    # --------------------------------------------------------------------------------------------------
                    # CONDENSATION
                    # --------------------------------------------------------------------------------------------------
                    m_cell_cell = element_stiffness_matrix[:_c0_c, :_c0_c]
                    m_cell_faces = element_stiffness_matrix[:_c0_c, _c0_c:]
                    m_faces_cell = element_stiffness_matrix[_c0_c:, :_c0_c]
                    m_faces_faces = element_stiffness_matrix[_c0_c:, _c0_c:]
                    # v_cell = element_residual[:_c0_c]
                    # v_faces = element_residual[_c0_c:]
                    v_cell = -element_residual[:_c0_c]
                    v_faces = -element_residual[_c0_c:]
                    m_cell_cell_inv = np.linalg.inv(m_cell_cell)
                    if debug_mode == 0:
                        print("TANGENT MATRIX IN CONDENSATION COND :")
                        print("{}".format(np.linalg.cond(m_cell_cell)))
                    # ge = m_faces_cell @ m_cell_cell_inv
                    # gd = (m_faces_cell @ m_cell_cell_inv) @ m_cell_faces
                    K_cond = m_faces_faces - ((m_faces_cell @ m_cell_cell_inv) @ m_cell_faces)
                    R_cond = v_faces - (m_faces_cell @ m_cell_cell_inv) @ v_cell
                    # --- SET CONDENSATION/DECONDENSATION MATRICES
                    # element.m_cell_cell_inv = m_cell_cell_inv
                    # element.m_cell_faces = m_cell_faces
                    # element.v_cell = v_cell
                    # --- ASSEMBLY
                    for _i_local, _i_global in enumerate(element.faces_indices):
                        _rg0 = _i_global * (_fk * _dx)
                        _rg1 = (_i_global + 1) * (_fk * _dx)
                        _re0 = _i_local * (_fk * _dx)
                        _re1 = (_i_local + 1) * (_fk * _dx)
                        residual[_rg0:_rg1] += R_cond[_re0:_re1]
                        for _j_local, _j_global in enumerate(element.faces_indices):
                            _cg0 = _j_global * (_fk * _dx)
                            _cg1 = (_j_global + 1) * (_fk * _dx)
                            _ce0 = _j_local * (_fk * _dx)
                            _ce1 = (_j_local + 1) * (_fk * _dx)
                            tangent_matrix[_rg0:_rg1, _cg0:_cg1] += K_cond[_re0:_re1, _ce0:_ce1]
                    # --- SET EXTERNAL FORCES COEFFICIENT
                    if np.max(np.abs(element_external_forces)) > external_forces_coefficient:
                        external_forces_coefficient = np.max(np.abs(element_external_forces))
                    if check and iteration > 0:
                        numerical_element_stiffness_matrix = (
                            numerical_forces_as_0[_element_index] - numerical_forces_as_1[_element_index]
                        ) / (2.0 * noise)
                        check_num = np.max(
                            np.abs(numerical_element_stiffness_matrix - element_stiffness_matrix)
                        ) / np.max(np.abs(element_stiffness_matrix))
                        if check_num > noise_tolerance:
                            print(
                                "AT ELEM : {} | CHECK_BEFORE_STAB : {}".format(str(_element_index).zfill(6), check_num)
                            )
                    if verbose:
                        print("ELEM : {} | INTERNAL_FORCES_AFTER STAB : \n {}".format(_element_index, element_internal_forces))
                        _iv0 = _element_index * len(element.cell.quadrature_weights)
                        _iv1 = (_element_index + 1) * len(element.cell.quadrature_weights)
                        print("ELEM : {} | DISPLACEMENT S0 : \n {}".format(_element_index, element.get_element_unknown_vector(self.field, self.finite_element, faces_unknown_vector_previous_step)))
                        print("ELEM : {} | GRADIENTS S0 : \n {}".format(_element_index, material.mat_data.s0.gradients[_iv0:_iv1]))
                        print("ELEM : {} | DISPLACEMENT S1 : \n {}".format(_element_index, element.get_element_unknown_vector(self.field, self.finite_element, faces_unknown_vector)))
                        print("ELEM : {} | GRADIENTS S1 : \n {}".format(_element_index, material.mat_data.s1.gradients[_iv0:_iv1]))
                        if not material.behaviour_name == "Elasticity":
                            print("ELEM : {} | INTERNAL_STATE_VARIABLES S0 : \n {}".format(_element_index, material.mat_data.s0.internal_state_variables[_iv0:_iv1]))
                            print("ELEM : {} | INTERNAL_STATE_VARIABLES S1 : \n {}".format(_element_index, material.mat_data.s1.internal_state_variables[_iv0:_iv1]))
                    # --- BOUNDARY CONDITIONS
                    for boundary_condition in self.boundary_conditions:
                        # --- DISPLACEMENT CONDITIONS
                        if boundary_condition.boundary_type == BoundaryType.DISPLACEMENT:
                            for f_local, f_global in enumerate(element.faces_indices):
                                if (
                                    f_global
                                    in self.mesh.faces_boundaries_connectivity[boundary_condition.boundary_name]
                                ):
                                    _l0 = _system_size + iter_face_constraint * _fk
                                    _l1 = _system_size + (iter_face_constraint + 1) * _fk
                                    _c0 = _cl * _dx + (f_local * _dx * _fk) + boundary_condition.direction * _fk
                                    _c1 = _cl * _dx + (f_local * _dx * _fk) + (boundary_condition.direction + 1) * _fk
                                    _r0 = f_global * _fk * _dx + _fk * boundary_condition.direction
                                    _r1 = f_global * _fk * _dx + _fk * (boundary_condition.direction + 1)
                                    # -------
                                    # face_displacement = element_unknown_increment[_r0:_r1]
                                    # face_displacement = np.copy(element_unknown_increment[_c0:_c1])
                                    # face_displacement = element_unknown_increment[_c0:_c1]
                                    face_lagrange = faces_unknown_vector[_l0:_l1]
                                    face_displacement = faces_unknown_vector[_r0:_r1]
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
                                    if debug_mode == 0:
                                        print("FACE MASS MATRIX IN DIRICHLET BOUND COND :")
                                        print("{}".format(np.linalg.cond(_m_psi_psi_face)))
                                    imposed_face_displacement = _m_psi_psi_face_inv @ _v_face_imposed_displacement
                                    face_displacement_difference = face_displacement - imposed_face_displacement
                                    # --- LAGRANGE INTERNAL FORCES PART
                                    # element_internal_forces[_c0:_c1] += (
                                    #     normalization_lagrange_coefficient * face_lagrange
                                    # )
                                    # residual[_r0:_r1] += (
                                    #     normalization_lagrange_coefficient * face_lagrange
                                    # )
                                    residual[_r0:_r1] -= (
                                        normalization_lagrange_coefficient * face_lagrange
                                    )
                                    # --- LAGRANGE MULTIPLIERS PART
                                    # residual[_l0:_l1] += (
                                    #     normalization_lagrange_coefficient * face_displacement_difference
                                    # )
                                    residual[_l0:_l1] -= (
                                        normalization_lagrange_coefficient * face_displacement_difference
                                    )
                                    # --- LAGRANGE MATRIX PART
                                    tangent_matrix[_l0:_l1, _r0:_r1] += normalization_lagrange_coefficient * np.eye(_fk, dtype=real)
                                    tangent_matrix[_r0:_r1, _l0:_l1] += normalization_lagrange_coefficient * np.eye(_fk, dtype=real)
                                    # --- SET EXTERNAL FORCES COEFFICIENT
                                    lagrange_external_forces = (
                                        normalization_lagrange_coefficient * imposed_face_displacement
                                    )
                                    if np.max(np.abs(lagrange_external_forces)) > external_forces_coefficient:
                                        external_forces_coefficient = np.max(np.abs(lagrange_external_forces))
                                    iter_face_constraint += 1
                        elif boundary_condition.boundary_type == BoundaryType.PRESSURE:
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
                                        _r0 = f_global * _fk * _dx + _fk * boundary_condition.direction
                                        _r1 = f_global * _fk * _dx + _fk * (boundary_condition.direction + 1)
                                        # element_external_forces[_c0:_c1] += vf
                                        # residual[_r0:_r1] -= vf
                                        residual[_r0:_r1] += vf
                                        if np.max(np.abs(vf)) > external_forces_coefficient:
                                            external_forces_coefficient = np.max(np.abs(vf))
                # --------------------------------------------------------------------------------------------------
                # RESIDUAL EVALUATION
                # --------------------------------------------------------------------------------------------------
                if external_forces_coefficient == 0.0:
                    external_forces_coefficient = 1.0
                residual_evaluation = np.max(np.abs(residual)) / external_forces_coefficient
                print("ITER : {} | RES_MAX : {}".format(str(iteration).zfill(4), residual_evaluation))
                # print("ITER : {} | =====================================================================".format(str(iteration).zfill(4)))
                # residual_values.append(residual_evaluation)
                # if residual_evaluation < self.tolerance:
                if residual_evaluation < self.tolerance or iteration == self.number_of_iterations - 1:
                    # ----------------------------------------------------------------------------------------------
                    # UPDATE INTERNAL VARIABLES
                    # ----------------------------------------------------------------------------------------------
                    mgis_bv.update(material.mat_data)
                    print("ITERATIONS : {}".format(iteration + 1))
                    self.write_vertex_res_files(file_suffix, faces_unknown_vector)
                    self.write_quadrature_points_res_files(file_suffix, material, faces_unknown_vector)
                    faces_unknown_vector_previous_step += faces_unknown_vector
                    residual_values.append(residual_evaluation)
                    break
                else:
                    # ----------------------------------------------------------------------------------------------
                    # SOLVE SYSTEM
                    # ----------------------------------------------------------------------------------------------
                    if debug_mode == 0:
                        print("K GLOBAL COND")
                        print("{}".format(np.linalg.cond(tangent_matrix)))
                    # sparse_global_matrix = csr_matrix(-tangent_matrix)
                    sparse_global_matrix = csr_matrix(tangent_matrix)
                    correction = spsolve(sparse_global_matrix, residual)
                    print("CORRECTION :")
                    # print(correction)
                    # correction = np.linalg.solve(-tangent_matrix, residual)
                    faces_unknown_vector += correction
                    if verbose:
                        print("R_K_RES : \n {}".format(residual - tangent_matrix @ correction))
                        print("R_K_RES_NORMALIZED : \n {}".format((residual - tangent_matrix @ correction)/external_forces_coefficient))
        plt.plot(range(len(residual_values)), residual_values, label="residual")
        plt.plot(range(len(residual_values)), [self.tolerance for local_i in range(len(residual_values))], label="tolerance")
        plt.ylabel("resiudal after convergence")
        plt.xlabel("time step index")
        plt.title("evolution of the residual")
        plt.legend()
        plt.gcf().subplots_adjust(left=0.15)
        plt.show()
        return