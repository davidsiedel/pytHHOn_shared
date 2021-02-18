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
    def test_square_traction_x_small_strain_linear_hardening(self):

        # ------------------------
        # POLYNOMIAL ORDER
        # ------------------------
        polynomial_order = 1

        # ------------------------
        # TIME
        # ------------------------
        spacing = 3
        time_steps = np.linspace(0.01, 0.05, spacing)
        # time_steps = [0.05]
        time_steps_1 = np.linspace(0.0, 7.e-3, spacing)
        time_steps_2 = np.linspace(7.e-3, -1.e-2, spacing)
        time_steps_3 = np.linspace(-1.e-2, 2.e-2, spacing)
        time_steps_4 = np.linspace(2.e-2, -3.e-2, spacing)
        time_steps_5 = np.linspace(-3.e-2, 4.e-2, spacing)
        time_steps = []
        for ts in [time_steps_1, time_steps_2[1:], time_steps_3[1:], time_steps_4[1:], time_steps_5[1:]]:
            # time_steps += list(np.sqrt(2.)*ts)
            time_steps += list(ts)
        print(time_steps)
        time_steps = np.array(time_steps)[:3]
        # time_steps = [0.0, 7.e-3, 1.e-2, 2.e-2, -3.e-2, 4.e-2]
        # time_steps = [0.0, 7.e-3, 1.e-2]

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
            # BoundaryCondition("TOP", pull, BoundaryType.PRESSURE, 1),
            # BoundaryCondition("LEFT", fixed, BoundaryType.DISPLACEMENT, 0),
            # BoundaryCondition("BOTTOM", fixed, BoundaryType.DISPLACEMENT, 1),
            # TRACTION
            BoundaryCondition("RIGHT", pull, BoundaryType.DISPLACEMENT, 0),
            BoundaryCondition("LEFT", fixed, BoundaryType.DISPLACEMENT, 0),
            BoundaryCondition("BOTTOM", fixed, BoundaryType.DISPLACEMENT, 1),
            # SHEAR
            # BoundaryCondition("TOP", pull, BoundaryType.DISPLACEMENT, 0),
            # BoundaryCondition("LEFT", fixed, BoundaryType.DISPLACEMENT, 1),
            # BoundaryCondition("BOTTOM", fixed, BoundaryType.DISPLACEMENT, 0),
            # BoundaryCondition("RIGHT", fixed, BoundaryType.DISPLACEMENT, 1),
        ]
        # ------------------------
        # MESH
        # ------------------------
        mesh_file_path = (
            "/home/dsiedel/Projects/pythhon/data/test_data/meshes/square_1.geof"
            # "/home/dsiedel/Projects/pythhon/data/test_data/meshes/square_5.geof"
            # "/home/dsiedel/Projects/pythhon/data/test_data/meshes/square_unsorted.geof"
            # "/home/dsiedel/Projects/pythhon/data/test_data/meshes/square_different.geof"
            # "/home/dsiedel/Projects/pythhon/data/test_data/meshes/square_different_2.geof"
            # "/home/dsiedel/Projects/pythhon/data/test_data/meshes/rectangle_test.geof"
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
        parameters = {"YoungModulus": 70.0e9, "PoissonRatio": 0.34, "HardeningSlope": 10.0e9, "YieldStress": 300.0e6}
        stabilization_parameter = parameters["YoungModulus"] / (1.0 + parameters["PoissonRatio"])
        # stabilization_parameter = 1.e4 * parameters["YoungModulus"]

        mat = Material(
            nq=p.mesh.number_of_cell_quadrature_points_in_mesh,
            library_path="/home/dsiedel/Projects/pythhon/behaviour/small_strain_isotropic_linear_hardening/src/libBehaviour.so",
            library_name="IsotropicLinearHardeningPlasticity",
            hypothesis=mgis_bv.Hypothesis.PLANESTRAIN,
            stabilization_parameter=stabilization_parameter,
            lagrange_parameter=parameters["YoungModulus"],
            field=displacement,
            parameters=None,
        )

        # --------------------------------------------------------------------------------------------------------------
        # LAUNCH
        # --------------------------------------------------------------------------------------------------------------
        reset_displacement_at_time_step = False
        # p.solve_newton_0(mat, reset_displacement_at_time_step)
        p.solve_newton_1(mat, reset_displacement_at_time_step)
        # p.solve_newton_check_1(mat, reset_displacement_at_time_step)
        # --------------------------------------------------------------------------------------------------------------
        # PROCESS RESULTS
        # --------------------------------------------------------------------------------------------------------------

        mtest_file_path = (
            "/home/dsiedel/Projects/pythhon/behaviour/testfront/traction_x_small_strain_isotropic_linear_hardening.res"
        )
        hho_file_path = "/home/dsiedel/Projects/pythhon/res"

        # plot_datat_2(mtest_file_path, hho_file_path, len(time_steps), 1, 5, 4, 8)

        plot_datat_2(mtest_file_path, hho_file_path, len(time_steps), 1, 5, 4, 8)
        plot_datat_2(mtest_file_path, hho_file_path, len(time_steps), 1, 6, 4, 9)
        plot_datat_2(mtest_file_path, hho_file_path, len(time_steps), 1, 7, 4, 10)
        plot_datat_2(mtest_file_path, hho_file_path, len(time_steps), 1, 8, 4, 11)

        # TRACTION COMPRESSION FINITE STRAIN
        # plot_data(mtest_file_path, hho_file_path, 1, 6, 1, 2, 11, 12)  # EPS_XX - SIG_XX
        # plot_data(mtest_file_path, hho_file_path, 1, 7, 1, 2, 13, 14)  # EPS_XX - SIG_YY
        # plot_data(mtest_file_path, hho_file_path, 1, 8, 1, 2, 15, 16)  # EPS_XX - SIG_ZZ
        # plot_data(mtest_file_path, hho_file_path, 1, 9, 1, 2, 17, 18)  # EPS_XX - SIG_XY

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

        # # ------------------------
        # # MATERIAL
        # # ------------------------
        # parameters = {
        #     "YoungModulus": 70.0e9,
        #     "PoissonRatio": 0.34,
        #     "HardeningSlope": 10.0e9,
        #     "YieldStress": 300.0e6,
        # }
        # stabilization_parameter = parameters["YoungModulus"] / (1. + parameters["PoissonRatio"])
        # # stabilization_parameter = parameters["YoungModulus"]
        #
        # mat = Material(
        #     p.mesh.number_of_cell_quadrature_points_in_mesh,
        #     "/home/dsiedel/Projects/pythhon/behaviour/small_strain_isotropic_linear_hardening/src/libBehaviour.so",
        #     "IsotropicLinearHardeningPlasticity",
        #     mgis_bv.Hypothesis.PLANESTRAIN,
        #     stabilization_parameter,
        #     field_type,
        #     parameters=None,
        #     # finite_strains=False
        # )

        # # ------------------------
        # # LAUNCH
        # # ------------------------
        # p.solve_newton_0(mat)
        #
        # from pp.plot_check_mtest import plot_data
        # mtest_file_path = "/home/dsiedel/Projects/pythhon/behaviour/testfront/small_strain_isotropic_linear_hardening.res"
        # hho_file_path = "/home/dsiedel/Projects/pythhon/res/res.csv"
        #
        # # TRACTION COMPRESSION
        # plot_data(mtest_file_path, hho_file_path, 1, 5, 1, 2, 9, 10)  # EPS_XX - SIG_XX
        # plot_data(mtest_file_path, hho_file_path, 1, 6, 1, 2, 11, 12)  # EPS_XX - SIG_YY
        # plot_data(mtest_file_path, hho_file_path, 1, 7, 1, 2, 13, 14)  # EPS_XX - SIG_ZZ
        # plot_data(mtest_file_path, hho_file_path, 1, 8, 1, 2, 15, 16)  # EPS_XX - SIG_XY

        self.assertTrue(True)
