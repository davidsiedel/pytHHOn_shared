from unittest import TestCase

import matplotlib
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FormatStrFormatter
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
    def test_ssna_finite_strain(self):

        # ------------------------
        # TIME
        # ------------------------

        time_steps = np.linspace(0.0, 6.e-3, 150)
        iterations = 20
        print(time_steps)

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
            return time

        def fixed(time: float, position: ndarray) -> float:
            return 0.0

        boundary_conditions = [
            # BoundaryCondition("LEFT", pull, BoundaryType.DISPLACEMENT, 1),
            # BoundaryCondition("RIGHT", fixed, BoundaryType.DISPLACEMENT, 0),
            # BoundaryCondition("RIGHT", fixed, BoundaryType.DISPLACEMENT, 1),
            # TRACTION
            # BoundaryCondition("RIGHT", pull, BoundaryType.DISPLACEMENT, 0),
            # BoundaryCondition("LEFT", fixed, BoundaryType.DISPLACEMENT, 0),
            # BoundaryCondition("BOTTOM", fixed, BoundaryType.DISPLACEMENT, 1),
            # TRACTION
            BoundaryCondition("LRU", pull, BoundaryType.DISPLACEMENT, 1),
            BoundaryCondition("LRD", fixed, BoundaryType.DISPLACEMENT, 1),
            BoundaryCondition("AXESYM", fixed, BoundaryType.DISPLACEMENT, 0),
            # TRACTION
            # BoundaryCondition("TOP", pull, BoundaryType.DISPLACEMENT, 1),
            # BoundaryCondition("BOTTOM", fixed, BoundaryType.DISPLACEMENT, 1),
            # BoundaryCondition("RIGHT", fixed, BoundaryType.DISPLACEMENT, 0),
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
            "meshes/ssna.geof"
        )
        # --------------------------------------------------------------------------------------------------------------
        # FIELD
        # --------------------------------------------------------------------------------------------------------------
        displacement = Field(
            label="U",
            euclidean_dimension=2,
            # field_type=FieldType.DISPLACEMENT_PLANE_STRAIN,
            field_type=FieldType.DISPLACEMENT_PLANE_STRESS,
            strain_type=StrainType.DISPLACEMENT_TRANSFORMATION_GRADIENT,
            stress_type=StressType.PIOLA_KIRCHOFF_1,
            derivation_type=DerivationType.FULL,
        )
        # displacement = Field(
        #     label="U",
        #     euclidean_dimension=2,
        #     field_type=FieldType.DISPLACEMENT_PLANE_STRAIN,
        #     strain_type=StrainType.DISPLACEMENT_SYMMETRIC_GRADIENT,
        #     stress_type=StressType.CAUCHY,
        #     derivation_type=DerivationType.SYMMETRIC,
        # )

        # --------------------------------------------------------------------------------------------------------------
        # FINITE ELEMENT
        # --------------------------------------------------------------------------------------------------------------
        finite_element = FiniteElement(
            element_type=ElementType.HDG_EQUAL,
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
            iterations=iterations,
            boundary_conditions=boundary_conditions,
            loads=loads,
            quadrature_type=QuadratureType.GAUSS,
            tolerance=1.e-4
        )

        # --------------------------------------------------------------------------------------------------------------
        # MATERIAL
        # --------------------------------------------------------------------------------------------------------------
        # parameters = {"YoungModulus": 70.e6, "PoissonRatio": 0.4999}
        # parameters = {"YoungModulus": 70.e6, "PoissonRatio": 0.4999, "HardeningSlope": 0.135e6, "YieldStress": 0.243e6}
        parameters = {"YoungModulus": 200.e9, "PoissonRatio": 0.3, "HardeningSlope": 20.e9, "YieldStress": 200.e6}
        parameters = {"YoungModulus": 70.e9, "PoissonRatio": 0.34, "HardeningSlope": 10.e9, "YieldStress": 300.e6}
        # parameters = {"YoungModulus": 70.0e9, "PoissonRatio": 0.34}
        # parameters = {"YoungModulus": 1.0, "PoissonRatio": 0.34}
        # parameters = {"YoungModulus": 1.0, "PoissonRatio": 0.0}
        # stabilization_parameter = parameters["YoungModulus"] / (1.0 + parameters["PoissonRatio"])
        stabilization_parameter = 1.0 * parameters["YoungModulus"] / (1.0 + parameters["PoissonRatio"])
        # stabilization_parameter = 70.e9 * parameters["YoungModulus"] / (1.0 + parameters["PoissonRatio"])
        # stabilization_parameter = 0. * parameters["YoungModulus"] / (1.0 + parameters["PoissonRatio"])
        # stabilization_parameter = 1. * parameters["YoungModulus"]

        mat = Material(
            nq=p.mesh.number_of_cell_quadrature_points_in_mesh,
            library_path="behaviour/src/libBehaviour.so",
            # library_name="IsotropicLinearHardeningPlasticity",
            # library_name="Elasticity",
            library_name="FiniteStrainIsotropicLinearHardeningPlasticity2",
            # hypothesis=mgis_bv.Hypothesis.PLANESTRAIN,
            hypothesis=mgis_bv.Hypothesis.PLANESTRESS,
            stabilization_parameter=stabilization_parameter,
            lagrange_parameter=parameters["YoungModulus"],
            # lagrange_parameter=1.0,
            field=displacement,
            parameters=None,
            # finite_strains=False
        )

        # --------------------------------------------------------------------------------------------------------------
        # LAUNCH
        # --------------------------------------------------------------------------------------------------------------
        # p.solve_newton_0(mat, reset_displacement_at_time_step)
        # p.solve_newton_1(mat)
        # p.solve_newton_1(mat)
        # p.solve_newton_2(mat, verbose=False, check=False)
        p.solve_newton_exact(mat, verbose=False, check=False)
        # p.solve_newton_check_1(mat, reset_displacement_at_time_step)
        # p.solve_newton_exact(mat, reset_displacement_at_time_step)
        # --------------------------------------------------------------------------------------------------------------
        # PROCESS RESULTS
        # --------------------------------------------------------------------------------------------------------------

        mtest_file_path = (
            "mtest/traction_isotropic_linear_hardening_finite_strain.res"
        )
        hho_file_path = "../../../res"

        res_folder = "../../../res"

        def __plot(column: int, time_step_index: int):
            from matplotlib.pyplot import figure
            figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
            _, _, filenames = next(walk(res_folder))
            for filename in filenames:
                if "{}".format(time_step_index).zfill(6) in filename and "qdp" in filename:
                    hho_file_path = path.join(res_folder, filename)
                    with open(hho_file_path, "r") as hho_res_file:
                        # fig, ax0d = plt.subplots(nrows=1, ncols=1)
                        fig, ax0d = plt.subplots(nrows=1, ncols=1, figsize=(11,15))
                        matplotlib.rcParams.update({'font.size': 22})
                        c_hho = hho_res_file.readlines()
                        field_label = c_hho[0].split(",")[column]
                        number_of_points = len(c_hho) - 1
                        eucli_d = displacement.euclidean_dimension
                        points = np.zeros((eucli_d, number_of_points), dtype=real)
                        field_vals = np.zeros((number_of_points,), dtype=real)
                        for l_count, line in enumerate(c_hho[1:]):
                            x_coordinates = float(line.split(",")[0])
                            y_coordinates = float(line.split(",")[1])
                            field_value = float(line.split(",")[column])
                            points[0, l_count] += x_coordinates
                            points[1, l_count] += y_coordinates
                            field_vals[l_count] += field_value
                        x, y = points
                        colors = [(0, 0, 1), (0, 1, 1), (0, 1, 0), (1, 1, 0), (1, 0, 0)]
                        perso = LinearSegmentedColormap.from_list("perso", colors, N=1000)
                        vmin = min(field_vals[:])
                        vmax = max(field_vals[:])
                        # levels = np.linspace(vmin, vmax, 50, endpoint=True)
                        levels = np.linspace(vmin, vmax, 1000, endpoint=True)
                        ticks = np.linspace(vmin, vmax, 10, endpoint=True)
                        datad = ax0d.tricontourf(x, y, field_vals[:], cmap=perso, levels=levels)
                        ax0d.get_xaxis().set_visible(False)
                        ax0d.get_yaxis().set_visible(False)
                        ax0d.set_xlabel("map of the domain $\Omega$")
                        cbar = fig.colorbar(datad, ax=ax0d, ticks=ticks)
                        cbar.set_label("{}".format(field_label), rotation=270, labelpad=15.0)
                        # plt.savefig("/home/dsiedel/Projects/pythhon/plots/{}.png".format(time_step))
                        plt.show()

        def __plot_det_f(time_step_index: int):
            _, _, filenames = next(walk(res_folder))
            for filename in filenames:
                if "{}".format(time_step_index).zfill(6) in filename and "qdp" in filename:
                    hho_file_path = path.join(res_folder, filename)
                    with open(hho_file_path, "r") as hho_res_file:
                        fig, ax0d = plt.subplots(nrows=1, ncols=1, figsize=(11,15))
                        matplotlib.rcParams.update({'font.size': 22})
                        c_hho = hho_res_file.readlines()
                        field_label = "DET_F"
                        number_of_points = len(c_hho) - 1
                        F = np.zeros((number_of_points, 3, 3), dtype=real)
                        eucli_d = displacement.euclidean_dimension
                        points = np.zeros((eucli_d, number_of_points), dtype=real)
                        field_vals = np.zeros((number_of_points,), dtype=real)
                        for l_count, line in enumerate(c_hho[1:]):
                            x_coordinates = float(line.split(",")[0])
                            y_coordinates = float(line.split(",")[1])
                            F_00 = float(line.split(",")[4])
                            F_11 = float(line.split(",")[5])
                            F_22 = float(line.split(",")[6])
                            F_01 = float(line.split(",")[7])
                            F_10 = float(line.split(",")[8])
                            F[l_count, 0, 0] = F_00
                            F[l_count, 1, 1] = F_11
                            F[l_count, 2, 2] = F_22
                            F[l_count, 0, 1] = F_01
                            F[l_count, 1, 0] = F_10
                            det_F = np.linalg.det(F[l_count])
                            points[0, l_count] += x_coordinates
                            points[1, l_count] += y_coordinates
                            field_vals[l_count] += det_F
                        x, y = points
                        colors = [(0, 0, 1), (0, 1, 1), (0, 1, 0), (1, 1, 0), (1, 0, 0)]
                        perso = LinearSegmentedColormap.from_list("perso", colors, N=1000)
                        vmin = min(field_vals[:])
                        vmax = max(field_vals[:])
                        # vmin = 9.84e-1
                        # vmax = 1.04
                        levels = np.linspace(vmin, vmax, 1000, endpoint=True)
                        ticks = np.linspace(vmin, vmax, 10, endpoint=True)
                        datad = ax0d.tricontourf(x, y, field_vals[:], cmap=perso, levels=levels)
                        ax0d.get_xaxis().set_visible(False)
                        ax0d.get_yaxis().set_visible(False)
                        ax0d.set_xlabel("map of the domain $\Omega$")
                        cbar = fig.colorbar(datad, ax=ax0d, ticks=ticks)
                        cbar.set_label("{}".format(field_label), rotation=270, labelpad=15.0)
                        # plt.savefig("/home/dsiedel/Projects/pythhon/plots/{}.png".format(time_step))
                        plt.show()

        def __plot_reaction_curve():
            _, _, filenames = next(walk(res_folder))
            forces = []
            for time_step_index in range(len(time_steps)):
                for filename in filenames:
                    if "{}".format(time_step_index).zfill(6) in filename and "qdp" in filename:
                        hho_file_path = path.join(res_folder, filename)
                        with open(hho_file_path, "r") as hho_res_file:
                            index = 10459
                            c_hho = hho_res_file.readlines()
                            line = c_hho[index]
                            x_coordinates = float(line.split(",")[0])
                            y_coordinates = float(line.split(",")[1])
                            x_disp = float(line.split(",")[2])
                            sig_11 = float(line.split(",")[10])
                            force = sig_11 * ((0.0054 + x_disp))
                            forces.append(force)
            cast_forces = []
            cast_times = []
            with open("castem/SSNA303_FU.csv", "r") as castfile:
                cast_c = castfile.readlines()
                for line in cast_c:
                    cast_force = float(line.split(",")[2])
                    cast_time = float(line.split(",")[1])
                    cast_forces.append(cast_force)
                    cast_times.append(cast_time)
            plt.plot(time_steps, forces, label="python HHO")
            plt.plot(cast_times, cast_forces, label="Cast3M", linestyle='--')
            plt.legend()
            plt.xlabel("displacement [m]")
            plt.ylabel("reaction force [N]")
            plt.grid()
            plt.show()
                            # fig, ax0d = plt.subplots(nrows=1, ncols=1)
                            # fig, ax0d = plt.subplots(nrows=1, ncols=1, figsize=(11,15))
                            # matplotlib.rcParams.update({'font.size': 22})
                            # c_hho = hho_res_file.readlines()
                            # field_label = c_hho[0].split(",")[column]
                            # number_of_points = len(c_hho) - 1
                            # eucli_d = displacement.euclidean_dimension
                            # points = np.zeros((eucli_d, number_of_points), dtype=real)
                            # field_vals = np.zeros((number_of_points,), dtype=real)
                            # for l_count, line in enumerate(c_hho[1:]):
                            #     x_coordinates = float(line.split(",")[0])
                            #     y_coordinates = float(line.split(",")[1])
                            #     field_value = float(line.split(",")[column])
                            #     points[0, l_count] += x_coordinates
                            #     points[1, l_count] += y_coordinates
                            #     field_vals[l_count] += field_value
                            # x, y = points
                            # colors = [(0, 0, 1), (0, 1, 1), (0, 1, 0), (1, 1, 0), (1, 0, 0)]
                            # perso = LinearSegmentedColormap.from_list("perso", colors, N=1000)
                            # vmin = min(field_vals[:])
                            # vmax = max(field_vals[:])
                            # # levels = np.linspace(vmin, vmax, 50, endpoint=True)
                            # levels = np.linspace(vmin, vmax, 1000, endpoint=True)
                            # ticks = np.linspace(vmin, vmax, 10, endpoint=True)
                            # datad = ax0d.tricontourf(x, y, field_vals[:], cmap=perso, levels=levels)
                            # ax0d.get_xaxis().set_visible(False)
                            # ax0d.get_yaxis().set_visible(False)
                            # ax0d.set_xlabel("map of the domain $\Omega$")
                            # cbar = fig.colorbar(datad, ax=ax0d, ticks=ticks)
                            # cbar.set_label("{}".format(field_label), rotation=270, labelpad=15.0)
                            # # plt.savefig("/home/dsiedel/Projects/pythhon/plots/{}.png".format(time_step))
                            # plt.show()

        # __plot(15, 149)
        # __plot(15, 149)
        __plot_det_f(149)
        # __plot_reaction_curve()

        # def plot_data_2(
        #         mtest_file_path: str,
        #         hho_res_dir_path: str,
        #         number_of_time_steps: int,
        #         m_x_inedx: int,
        #         m_y_index: int,
        #         d_x_inedx: int,
        #         d_y_inedx: int,
        # ):
        #     coef = 1.0e6
        #     with open(mtest_file_path, "r") as mres_file:
        #         c = mres_file.readlines()
        #         labels = []
        #         m_x_data = []
        #         m_y_data = []
        #         for line in c:
        #             if line[0] == "#":
        #                 alpha = line.split("column:")
        #                 labels.append(alpha[1].replace("\n", ""))
        #         line_start = len(labels)
        #         for time, line in enumerate(c[line_start:]):
        #             m_x_data.append(float(line.split(" ")[m_x_inedx]))
        #             m_y_data.append(float(line.split(" ")[m_y_index]) / coef)
        #         m_x_label = labels[m_x_inedx]
        #         m_y_label = labels[m_y_index]
        #     _, _, filenames = next(walk(hho_res_dir_path))
        #     d_x_min_data = []
        #     d_y_min_data = []
        #     d_x_max_data = []
        #     d_y_max_data = []
        #     for time_step_index in range(number_of_time_steps):
        #         for filename in filenames:
        #             if "{}".format(time_step_index).zfill(6) in filename and "qdp" in filename:
        #                 hho_file_path = path.join(hho_res_dir_path, filename)
        #                 with open(hho_file_path, "r") as hho_res_file:
        #                     c_hho = hho_res_file.readlines()
        #                     x_val_min = +np.inf
        #                     x_val_max = -np.inf
        #                     y_val_min = +np.inf
        #                     y_val_max = -np.inf
        #                     for line in c_hho[1:]:
        #                         x_val = float(line.split(",")[d_x_inedx])
        #                         y_val = float(line.split(",")[d_y_inedx])
        #                         if x_val > x_val_max:
        #                             x_val_max = x_val
        #                         if x_val < x_val_min:
        #                             x_val_min = x_val
        #                         if y_val > y_val_max:
        #                             y_val_max = y_val
        #                         if y_val < y_val_min:
        #                             y_val_min = y_val
        #                     d_x_min_data.append(x_val_min)
        #                     d_y_min_data.append(y_val_min / coef)
        #                     d_x_max_data.append(x_val_max)
        #                     d_y_max_data.append(y_val_max / coef)
        #     # fig, ax = plt.subplots()
        #     plt.plot(m_x_data, m_y_data, color="blue", label="MTEST")
        #     plt.plot(d_x_min_data, d_y_min_data, color="green", linestyle="--", label="HHO_MIN")
        #     plt.plot(d_x_max_data, d_y_max_data, color="purple", linestyle="--", label="HHO_MAX")
        #     plt.xlabel(m_x_label)
        #     plt.ylabel(m_y_label + " [MPA]")
        #     # plt.title("MTEST HHO COMPARISON CYCLIC LOADING ISOTROPIC LINEAR HARDENING")
        #     plt.legend()
        #     plt.gcf().subplots_adjust(left=0.15)
        #     plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        #     # ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        #     # plt.ylim(-0.35, 0.35) #TRACTION
        #     # plt.ylim(-0.01* 1.e-8, 0.009 * 1.e-7)
        #     # plt.ylim(-50.0, 1000.0)  # SHEAR
        #     # plt.xlim(-0.05, 0.05)
        #     plt.grid()
        #     plt.show()
        #     return
        #
        # # CONTRAINTES
        # # plot_data_2(mtest_file_path, hho_file_path, len(time_steps), 1, 5, 4, 8)
        # # plot_data_2(mtest_file_path, hho_file_path, len(time_steps), 1, 6, 4, 9)
        # # plot_data_2(mtest_file_path, hho_file_path, len(time_steps), 1, 7, 4, 10)
        # # plot_data_2(mtest_file_path, hho_file_path, len(time_steps), 1, 8, 4, 11)
        #
        # plot_data_2(mtest_file_path, hho_file_path, len(time_steps), 1, 6, 4, 9)
        # plot_data_2(mtest_file_path, hho_file_path, len(time_steps), 1, 7, 4, 10)
        # plot_data_2(mtest_file_path, hho_file_path, len(time_steps), 1, 8, 4, 11)
        # plot_data_2(mtest_file_path, hho_file_path, len(time_steps), 1, 9, 4, 12)
        # #
        # # plot_data_2(mtest_file_path, hho_file_path, len(time_steps), 1, 1, 4, 4)
        # # plot_data_2(mtest_file_path, hho_file_path, len(time_steps), 1, 2, 4, 5)
        # # plot_data_2(mtest_file_path, hho_file_path, len(time_steps), 1, 3, 4, 6)
        # # plot_data_2(mtest_file_path, hho_file_path, len(time_steps), 1, 4, 4, 7)

        self.assertTrue(True)
