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
    def test_cook_small_strain_linear_hardening(self):
        # ------------------------
        # POLYNOMIAL ORDER
        # ------------------------
        polynomial_order = 1

        # ------------------------
        # TIME
        # ------------------------
        # F_min = 0.000
        # F_max = 1./16.
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
        ]
        # ------------------------
        # MESH
        # ------------------------
        mesh_file_path = (
            # "/home/dsiedel/Projects/pythhon/data/test_data/meshes/cook_1.geof"
            "/home/dsiedel/Projects/pythhon/data/test_data/meshes/cook5.geof"
            # "/home/dsiedel/Projects/pythhon/data/test_data/meshes/cook_20.geof"
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
        # stabilization_parameter = 1.e3 * parameters["YoungModulus"]

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
        # mat_0 = Material(
        #     nq=p.mesh.number_of_cell_quadrature_points_in_mesh,
        #     library_path="/home/dsiedel/Projects/pythhon/behaviour/small_strain_isotropic_linear_hardening/src/libBehaviour.so",
        #     library_name="IsotropicLinearHardeningPlasticity",
        #     hypothesis=mgis_bv.Hypothesis.PLANESTRAIN,
        #     stabilization_parameter=stabilization_parameter,
        #     lagrange_parameter=parameters["YoungModulus"],
        #     field=displacement,
        #     parameters=None,
        # )
        # mat_1 = Material(
        #     nq=p.mesh.number_of_cell_quadrature_points_in_mesh,
        #     library_path="/home/dsiedel/Projects/pythhon/behaviour/small_strain_isotropic_linear_hardening/src/libBehaviour.so",
        #     library_name="IsotropicLinearHardeningPlasticity",
        #     hypothesis=mgis_bv.Hypothesis.PLANESTRAIN,
        #     stabilization_parameter=stabilization_parameter,
        #     lagrange_parameter=parameters["YoungModulus"],
        #     field=displacement,
        #     parameters=None,
        # )

        # --------------------------------------------------------------------------------------------------------------
        # LAUNCH
        # --------------------------------------------------------------------------------------------------------------
        reset_displacement_at_time_step = False
        # p.solve_newton_0(mat, reset_displacement_at_time_step)
        # p.solve_newton_check(mat, mat_0, mat_1, reset_displacement_at_time_step)
        # p.solve_newton_1(mat, reset_displacement_at_time_step)
        p.solve_newton_check_1(mat, reset_displacement_at_time_step)
        # p.solve_newton_2(mat, reset_displacement_at_time_step)

        self.assertTrue(True)
