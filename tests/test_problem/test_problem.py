from unittest import TestCase

from numpy import ndarray

from pythhon.pbbb.problem import Problem
from pythhon.pbbb.boundary_condition import BoundaryCondition
from pythhon.pbbb.load import Load
from pythhon.parameters import *
from data.test_data.meshes.config import test_parsers_folder_path

from mgis import behaviour as mgis_bv

from pythhon.pbbb.material import Material

class TestProblem(TestCase):
    def test_build(self):
        # ------------------------
        # MESH
        # ------------------------
        # mesh_file_path = os.path.join(test_parsers_folder_path, "c2d3_2.geof")
        mesh_file_path = os.path.join(test_parsers_folder_path, "cook5.geof")
        # mesh_file_path = os.path.join(test_parsers_folder_path, "cook80.geof")
        # mesh_file_path = os.path.join(test_parsers_folder_path, "cook40.geof")

        # ------------------------
        # POLYNOMIAL ORDER
        # ------------------------
        polynomial_order = 1

        # ------------------------
        # TIME
        # ------------------------
        time_steps = np.linspace(1.0 / 32.0, 1.0 / 16.0, 8)
        iterations = 200

        # ------------------------
        # LOAD
        # ------------------------
        def load_x(time: float, position: ndarray) -> float:
            """

            :param time:
            :param position: in d -dimension
            :return:
            """
            f_x = -2.0 * (np.pi ** 2) * (np.sin(np.pi * position[0])) * (np.sin(np.pi * position[1]))
            return f_x

        # ------------------------
        # BC
        # ------------------------
        def pull(time: float, position: ndarray) -> float:
            """

            :param time:
            :param position: in d-1 dimension
            :return:
            """
            return time

        def fixed(time: float, position: ndarray) -> float:
            """

            :param time:
            :param position: in d-1 dimension
            :return:
            """
            return 0.0

        boundary_conditions = [
            BoundaryCondition("RIGHT", pull, BoundaryType.PRESSURE, 1),
            # BoundaryCondition("RIGHT", pull, BoundaryType.DISPLACEMENT, 1),
            # BoundaryCondition("RIGHT", pull, BoundaryType.DISPLACEMENT, 0),
            BoundaryCondition("LEFT", fixed, BoundaryType.DISPLACEMENT, 0),
            BoundaryCondition("LEFT", fixed, BoundaryType.DISPLACEMENT, 1),
        ]
        # boundary_conditions = [
        #     BoundaryCondition("TOP", pull, BoundaryType.DISPLACEMENT, 1),
        #     BoundaryCondition("BOTTOM", fixed, BoundaryType.DISPLACEMENT, 0),
        # ]
        loads = [
            Load(load_x, 0),
            # Load(load_x, 1),
        ]
        loads = None

        # ------------------------
        # PROBLEM DATA
        # ------------------------
        quadrature_type: QuadratureType = QuadratureType.GAUSS
        basis_type: BasisType = BasisType.MONOMIAL
        element_type: ElementType = ElementType.HDG_EQUAL
        field_type: FieldType = FieldType.DISPLACEMENT_SMALL_STRAIN_PLANE_STRAIN

        # ------------------------
        # LAUCH
        # ------------------------
        p = Problem(
            mesh_file_path,
            polynomial_order,
            time_steps,
            iterations,
            boundary_conditions,
            loads=loads,
            quadrature_type=quadrature_type,
            basis_type=basis_type,
            element_type=element_type,
            field_type=field_type,
            # gradient_type=gradient_type,
        )

        material_params = {
            "YoungModulus": 1.124999981250001,
            "PoissonRatio": 0.49999999,
        }
        external_params = {
            "Temperature": 293.15,
        }

        material = Material(
            "/behaviour/trash/src/libBehaviour.so",
            "Elasticity",
            mgis_bv.Hypothesis.PLANESTRAIN,
            material_params,
            external_params,
        )

        p.solve_newton(
            material
        )
        self.assertTrue(True)
