from unittest import TestCase

from numpy import ndarray

from pythhon.pbbb.problem import Problem
from pythhon.pbbb.boundary_condition import BoundaryCondition
from pythhon.pbbb.load import Load
from pythhon.pbbb.field import Field
from pythhon.fem.element.finite_element import FiniteElement
from pythhon.parameters import *
from pp.post_processing import *

from mgis import behaviour as mgis_bv

from pythhon.pbbb.material import Material


class TestMecha(TestCase):
    def test_cook_elasticity(self):
        # ------------------------
        # POLYNOMIAL ORDER
        # ------------------------
        polynomial_order = 1

        # ------------------------
        # TIME
        # ------------------------
        # F_min = 0.001
        # F_max = 1./16.
        # F_min = 0.0 / 16.
        # F_max = 70.e9 / 16.
        # time_steps = np.linspace(F_min, F_max, 20)
        # iterations = 20
        P_min = 0.0
        P_max = 70.e9 / 16.
        time_steps = np.linspace(P_min, P_max, 100)
        iterations = 100

        # ------------------------
        # LOAD
        # ------------------------
        def volumetric_load(time: float, position: ndarray):
            return 0

        loads = [Load(volumetric_load, 0), Load(volumetric_load, 1)]

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
            BoundaryCondition("LEFT", fixed, BoundaryType.DISPLACEMENT, 0),
            BoundaryCondition("LEFT", fixed, BoundaryType.DISPLACEMENT, 1),
            #
            # BoundaryCondition("RIGHT", pull, BoundaryType.DISPLACEMENT, 0),
            # BoundaryCondition("LEFT", fixed, BoundaryType.DISPLACEMENT, 0),
            # BoundaryCondition("BOTTOM", fixed, BoundaryType.DISPLACEMENT, 1),
        ]
        # ------------------------
        # MESH
        # ------------------------
        mesh_file_path = (
            # "/home/dsiedel/Projects/pythhon/data/test_data/meshes/cook_1.geof"
            # "/home/dsiedel/Projects/pythhon/data/test_data/meshes/cook5.geof"
            "/home/dsiedel/Projects/pythhon/data/test_data/meshes/cook_20.geof"
        )
        # --------------------------------------------------------------------------------------------------------------
        # FIELD
        # --------------------------------------------------------------------------------------------------------------
        displacement = Field(
            label="U",
            euclidean_dimension=2,
            field_type=FieldType.DISPLACEMENT_PLANE_STRAIN,
            strain_type=StrainType.DISPLACEMENT_SYMMETRIC_GRADIENT,
            stress_type=StressType.CAUCHY,
            derivation_type=DerivationType.SYMMETRIC,
        )

        # --------------------------------------------------------------------------------------------------------------
        # FINITE ELEMENT
        # --------------------------------------------------------------------------------------------------------------
        finite_element = FiniteElement(
            element_type=ElementType.HDG_EQUAL,
            polynomial_order=polynomial_order,
            euclidean_dimension=2,
            basis_type=BasisType.MONOMIAL
        )

        # --------------------------------------------------------------------------------------------------------------
        # PROBLEM
        # --------------------------------------------------------------------------------------------------------------
        p = Problem(
            mesh_file_path=mesh_file_path,
            field=displacement,
            polynomial_order=polynomial_order,
            finite_element=finite_element,
            time_steps=time_steps,
            iterations=iterations,
            boundary_conditions=boundary_conditions,
            loads=loads,
            quadrature_type=QuadratureType.GAUSS,
            tolerance=1.0e-4
        )

        # --------------------------------------------------------------------------------------------------------------
        # MATERIAL
        # --------------------------------------------------------------------------------------------------------------
        parameters = {"YoungModulus": 70.0e9, "PoissonRatio": 0.34}
        stabilization_parameter = parameters["YoungModulus"] / (1.0 + parameters["PoissonRatio"])
        # stabilization_parameter = 1.e2 * parameters["YoungModulus"]

        mat = Material(
            nq=p.mesh.number_of_cell_quadrature_points_in_mesh,
            library_path="/home/dsiedel/Projects/pythhon/behaviour/elasticity/src/libBehaviour.so",
            library_name="Elasticity",
            hypothesis=mgis_bv.Hypothesis.PLANESTRAIN,
            stabilization_parameter=stabilization_parameter,
            lagrange_parameter=parameters["YoungModulus"],
            field=displacement,
            parameters=None,
            # finite_strains=False
        )

        # --------------------------------------------------------------------------------------------------------------
        # LAUNCH
        # --------------------------------------------------------------------------------------------------------------
        reset_displacement_at_time_step = True
        # p.solve_newton_0(mat, reset_displacement_at_time_step)
        # p.solve_newton_2(mat, reset_displacement_at_time_step)
        p.solve_newton_check_1(mat, reset_displacement_at_time_step)

        # --------------------------------------------------------------------------------------------------------------
        # PROCESS RESULTS
        # --------------------------------------------------------------------------------------------------------------

        # mtest_file_path = (
        #     "/home/dsiedel/Projects/pythhon/behaviour/testfront/elasticity.res"
        # )
        # hho_file_path = "/home/dsiedel/Projects/pythhon/res"
        #
        # plot_datat_2(mtest_file_path, hho_file_path, len(time_steps), 1, 5, 4, 8)


        # # ------------------------
        # # MESH
        # # ------------------------
        # mesh_file_path = (
        #     "/home/dsiedel/Projects/pythhon/data/test_data/meshes/cook5.geof"
        # )
        #
        # # ------------------------
        # # POLYNOMIAL ORDER
        # # ------------------------
        # polynomial_order = 1
        #
        # # ------------------------
        # # TIME
        # # ------------------------
        # # F_min = 0.000
        # # F_max = 1./16.
        # F_min = 3.e9/16.
        # F_max = 70.e9/16.
        # time_steps = np.linspace(F_min, F_max, 8)
        # iterations = 20
        #
        # # ------------------------
        # # LOAD
        # # ------------------------
        # loads = None
        #
        # # ------------------------
        # # BC
        # # ------------------------
        # def pull(time: float, position: ndarray) -> float:
        #     """
        #
        #     :param time:
        #     :param position: in d-1 dimension
        #     :return:
        #     """
        #     return time
        #
        # def fixed(time: float, position: ndarray) -> float:
        #     """
        #
        #     :param time:
        #     :param position: in d-1 dimension
        #     :return:
        #     """
        #     return 0.0
        #
        # boundary_conditions = [
        #     BoundaryCondition("RIGHT", pull, BoundaryType.PRESSURE, 1),
        #     BoundaryCondition("LEFT", fixed, BoundaryType.DISPLACEMENT, 0),
        #     BoundaryCondition("LEFT", fixed, BoundaryType.DISPLACEMENT, 1),
        # ]
        #
        # # ------------------------
        # # PROBLEM DATA
        # # ------------------------
        # quadrature_type: QuadratureType = QuadratureType.GAUSS
        # basis_type: BasisType = BasisType.MONOMIAL
        # element_type: ElementType = ElementType.HDG_EQUAL
        # field_type: FieldType = FieldType.DISPLACEMENT_SMALL_STRAIN_PLANE_STRAIN
        #
        # # ------------------------
        # # DEFINE PROBLEM
        # # ------------------------
        # p = Problem(
        #     mesh_file_path,
        #     polynomial_order,
        #     time_steps,
        #     iterations,
        #     boundary_conditions,
        #     loads=loads,
        #     quadrature_type=quadrature_type,
        #     basis_type=basis_type,
        #     element_type=element_type,
        #     field_type=field_type,
        #     tolerance=1.e-4
        # )
        #
        # # ------------------------
        # # MATERIAL
        # # ------------------------
        # parameters = {
        #     # "YoungModulus": 1.0,
        #     # "PoissonRatio": 0.4999,
        #     "YoungModulus": 70.0e9,
        #     "PoissonRatio": 0.34,
        # }
        # stabilization_parameter = parameters["YoungModulus"] / (1. + parameters["PoissonRatio"])
        #
        # mat = Material(
        #     p.mesh.number_of_cell_quadrature_points_in_mesh,
        #     "/home/dsiedel/Projects/pythhon/behaviour/elasticity/src/libBehaviour.so",
        #     "Elasticity",
        #     mgis_bv.Hypothesis.PLANESTRAIN,
        #     stabilization_parameter,
        #     field_type,
        #     parameters=parameters,
        #     # finite_strains=False
        #
        # )
        #
        # # ------------------------
        # # LAUNCH
        # # ------------------------
        # reset_displacement_at_time_step = True
        # p.solve_newton_0(mat, reset_displacement_at_time_step)

        self.assertTrue(True)
