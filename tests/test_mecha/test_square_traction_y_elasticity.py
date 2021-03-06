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
    def test_square_traction_y_elsaticity(self):
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
        time_steps = np.array(time_steps)
        # time_steps = [0.0, 7.e-3, -1.e-2, 2.e-2, -3.e-2, 4.e-2]

        iterations = 20

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
            # BoundaryCondition("LEFT", pull, BoundaryType.PRESSURE, 1),
            # BoundaryCondition("RIGHT", fixed, BoundaryType.DISPLACEMENT, 0),
            # BoundaryCondition("RIGHT", fixed, BoundaryType.DISPLACEMENT, 1),
            # TRACTION
            # BoundaryCondition("RIGHT", pull, BoundaryType.DISPLACEMENT, 0),
            # BoundaryCondition("LEFT", fixed, BoundaryType.DISPLACEMENT, 0),
            # BoundaryCondition("BOTTOM", fixed, BoundaryType.DISPLACEMENT, 1),
            # TRACTION
            BoundaryCondition("TOP", pull, BoundaryType.DISPLACEMENT, 1),
            BoundaryCondition("BOTTOM", fixed, BoundaryType.DISPLACEMENT, 1),
            BoundaryCondition("RIGHT", fixed, BoundaryType.DISPLACEMENT, 0),
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
            # "/home/dsiedel/Projects/pythhon/data/test_data/meshes/square_1.geof"
            # "/home/dsiedel/Projects/pythhon/data/test_data/meshes/square_5.geof"
            "/home/dsiedel/Projects/pythhon/data/test_data/meshes/square_unsorted.geof"
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
            # field=field,
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
        # stabilization_parameter = 1. * parameters["YoungModulus"]

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
        # p.solve_newton_1(mat, reset_displacement_at_time_step)
        p.solve_newton_check_1(mat, reset_displacement_at_time_step)
        # --------------------------------------------------------------------------------------------------------------
        # PROCESS RESULTS
        # --------------------------------------------------------------------------------------------------------------

        mtest_file_path = (
            "/home/dsiedel/Projects/pythhon/behaviour/testfront/traction_y_elasticity.res"
        )
        hho_file_path = "/home/dsiedel/Projects/pythhon/res"

        plot_datat_2(mtest_file_path, hho_file_path, len(time_steps), 2, 5, 5, 8)
        plot_datat_2(mtest_file_path, hho_file_path, len(time_steps), 2, 6, 5, 9)
        plot_datat_2(mtest_file_path, hho_file_path, len(time_steps), 2, 7, 5, 10)
        plot_datat_2(mtest_file_path, hho_file_path, len(time_steps), 2, 8, 5, 11)

        self.assertTrue(True)
