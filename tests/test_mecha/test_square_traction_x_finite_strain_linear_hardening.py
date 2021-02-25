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
    def test_square_traction_x_finite_strain_linear_hardening(self):

        # ------------------------
        # POLYNOMIAL ORDER
        # ------------------------

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

        iterations = 7

        # time_steps = np.linspace(0.0, 0.008, 9)
        # time_steps = np.linspace(0.0, 0.004, 9)

        number_of_iterations = 10


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
            "../../data/test_data/meshes/square_1.geof"
            # "../../data/test_data/meshes/square_5.geof"
            # "../../data/test_data/meshes/square_unsorted.geof"
            # "../../data/test_data/meshes/square_different.geof"
            # "../../data/test_data/meshes/square_different_2.geof"
            # "../../data/test_data/meshes/rectangle_test.geof"
        )
        # --------------------------------------------------------------------------------------------------------------
        # FIELD
        # --------------------------------------------------------------------------------------------------------------
        displacement = Field(
            label="U",
            euclidean_dimension=2,
            field_type=FieldType.DISPLACEMENT_PLANE_STRAIN,
            strain_type=StrainType.DISPLACEMENT_TRANSFORMATION_GRADIENT,
            stress_type=StressType.PIOLA_KIRCHOFF_1,
            derivation_type=DerivationType.FULL,
        )

        # --------------------------------------------------------------------------------------------------------------
        # FINITE ELEMENT
        # --------------------------------------------------------------------------------------------------------------
        finite_element = FiniteElement(
            element_type=ElementType.HDG_EQUAL,
            # element_type=ElementType.HDG_HIGH,
            # element_type=ElementType.HDG_LOW,
            polynomial_order=1,
            euclidean_dimension=2,
            basis_type=BasisType.MONOMIAL
        )

        # --------------------------------------------------------------------------------------------------------------
        # PROBLEM
        # --------------------------------------------------------------------------------------------------------------
        p = Problem(
            mesh_file_path=mesh_file_path,
            field=displacement,
            polynomial_order=finite_element.polynomial_order,
            finite_element=finite_element,
            time_steps=time_steps,
            # iterations=iterations,
            iterations=number_of_iterations,
            boundary_conditions=boundary_conditions,
            loads=loads,
            quadrature_type=QuadratureType.GAUSS,
            tolerance=1.0e-6
        )

        # --------------------------------------------------------------------------------------------------------------
        # MATERIAL
        # --------------------------------------------------------------------------------------------------------------
        parameters = {"YoungModulus": 70.0e9, "PoissonRatio": 0.34, "HardeningSlope": 10.0e9, "YieldStress": 300.0e6}
        stabilization_parameter = parameters["YoungModulus"] / (1.0 + parameters["PoissonRatio"])
        # stabilization_parameter = 70e9 * parameters["YoungModulus"] / (1.0 + parameters["PoissonRatio"])
        # stabilization_parameter = 1.e2 * parameters["YoungModulus"]

        mat = Material(
            nq=p.mesh.number_of_cell_quadrature_points_in_mesh,
            library_path="../../behaviour/finite_strain_isotropic_linear_hardening/src/libBehaviour.so",
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
        # p.solve_newton_0(mat, reset_displacement_at_time_step)
        # p.solve_newton_1(mat, reset_displacement_at_time_step)
        # p.solve_newton_check_1(mat, reset_displacement_at_time_step)
        # p.solve_newton_2(mat, verbose=False, check=False)
        p.solve_newton_exact(mat, verbose=False, check=False)
        # --------------------------------------------------------------------------------------------------------------
        # PROCESS RESULTS
        # --------------------------------------------------------------------------------------------------------------

        mtest_file_path = (
            "../../behaviour/testfront/traction_x_finite_strain_isotropic_linear_hardening.res"
            # "../../behaviour/testfront/traction_x_finite_strain_isotropic_linear_hardening_2.res"
        )
        hho_file_path = "../../res"

        def plot_data(
                mtest_file_path: str,
                hho_res_dir_path: str,
                number_of_time_steps: int,
                m_x_inedx: int,
                m_y_index: int,
                d_x_inedx: int,
                d_y_inedx: int,
        ):
            coef = 10.0e6
            with open(mtest_file_path, "r") as mres_file:
                c = mres_file.readlines()
                labels = []
                m_x_data = []
                m_y_data = []
                for line in c:
                    if line[0] == "#":
                        alpha = line.split("column:")
                        labels.append(alpha[1].replace("\n", ""))
                line_start = len(labels)
                for time, line in enumerate(c[line_start:]):
                    m_x_data.append(float(line.split(" ")[m_x_inedx]))
                    m_y_data.append(float(line.split(" ")[m_y_index]) / coef)
                m_x_label = labels[m_x_inedx]
                m_y_label = labels[m_y_index]
            _, _, filenames = next(walk(hho_res_dir_path))
            d_x_min_data = []
            d_y_min_data = []
            d_x_max_data = []
            d_y_max_data = []
            for time_step_index in range(number_of_time_steps):
                for filename in filenames:
                    if "{}".format(time_step_index).zfill(6) in filename and "qdp" in filename:
                        hho_file_path = path.join(hho_res_dir_path, filename)
                        with open(hho_file_path, "r") as hho_res_file:
                            c_hho = hho_res_file.readlines()
                            x_val_min = +np.inf
                            x_val_max = -np.inf
                            y_val_min = +np.inf
                            y_val_max = -np.inf
                            for line in c_hho[1:]:
                                x_val = float(line.split(",")[d_x_inedx])
                                y_val = float(line.split(",")[d_y_inedx])
                                if x_val > x_val_max:
                                    x_val_max = x_val
                                if x_val < x_val_min:
                                    x_val_min = x_val
                                if y_val > y_val_max:
                                    y_val_max = y_val
                                if y_val < y_val_min:
                                    y_val_min = y_val
                            d_x_min_data.append(x_val_min)
                            d_y_min_data.append(y_val_min / coef)
                            d_x_max_data.append(x_val_max)
                            d_y_max_data.append(y_val_max / coef)
            plt.plot(m_x_data, m_y_data, color="blue", label="MTEST")
            plt.plot(d_x_min_data, d_y_min_data, color="green", linestyle="--", label="HHO_MIN")
            plt.plot(d_x_max_data, d_y_max_data, color="purple", linestyle="--", label="HHO_MAX")
            plt.xlabel(m_x_label)
            plt.ylabel(m_y_label + " [MPA]")
            # plt.title("MTEST HHO COMPARISON CYCLIC LOADING ISOTROPIC LINEAR HARDENING")
            plt.legend()
            # plt.ylim(-0.35, 0.35) #TRACTION
            # plt.ylim(-200.0, 200.0)  # SHEAR
            # plt.ylim(-400.0, 400.0)  # SHEAR
            # plt.xlim(-0.05, 0.05)
            plt.gcf().subplots_adjust(left=0.15)
            plt.grid()
            plt.show()
            return

        # plot_datat_2(mtest_file_path, hho_file_path, len(time_steps), 1, 5, 4, 8)

        plot_data(mtest_file_path, hho_file_path, len(time_steps), 1, 6, 4, 9)
        plot_data(mtest_file_path, hho_file_path, len(time_steps), 1, 7, 4, 10)
        plot_data(mtest_file_path, hho_file_path, len(time_steps), 1, 8, 4, 11)
        plot_data(mtest_file_path, hho_file_path, len(time_steps), 1, 9, 4, 12)
        # plot_datat_2(mtest_file_path, hho_file_path, len(time_steps), 1, 8, 4, 13)

        self.assertTrue(True)
