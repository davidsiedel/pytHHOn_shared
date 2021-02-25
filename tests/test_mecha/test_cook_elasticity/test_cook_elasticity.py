from unittest import TestCase

from matplotlib.colors import LinearSegmentedColormap
from numpy import ndarray

from pythhon.pbbb.problem import Problem
from pythhon.pbbb.boundary_condition import BoundaryCondition
from pythhon.pbbb.load import Load
from pythhon.pbbb.field import Field
from pythhon.fem.element.finite_element import FiniteElement
from pythhon.parameters import *
from pp.post_processing import *

import matplotlib.pyplot as plt

from mgis import behaviour as mgis_bv

from pythhon.pbbb.material import Material


class TestMecha(TestCase):
    def test_cook_elasticity(self):
        # ------------------------
        # POLYNOMIAL ORDER
        # ------------------------

        # ------------------------
        # TIME
        # ------------------------
        P_min = 0.0
        P_max = 70.e9 / 16.
        P_min = 0.01
        P_max = 1. / 16.
        time_steps = np.linspace(P_min, P_max, 3)
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
            "meshes/cook_20.geof"
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
            tolerance=1.0e-6
        )

        # --------------------------------------------------------------------------------------------------------------
        # MATERIAL
        # --------------------------------------------------------------------------------------------------------------
        parameters = {"YoungModulus": 1.13, "PoissonRatio": 0.4999}
        stabilization_parameter = 1.e-19 * parameters["YoungModulus"] / (1.0 + parameters["PoissonRatio"])
        stabilization_parameter = 1.0 * parameters["YoungModulus"] / (1.0 + parameters["PoissonRatio"])
        # stabilization_parameter = 1.e2 * parameters["YoungModulus"]

        mat = Material(
            nq=p.mesh.number_of_cell_quadrature_points_in_mesh,
            library_path="behaviour/src/libBehaviour.so",
            library_name="Elasticity",
            hypothesis=mgis_bv.Hypothesis.PLANESTRAIN,
            stabilization_parameter=stabilization_parameter,
            lagrange_parameter=parameters["YoungModulus"],
            field=displacement,
            parameters=None,
        )

        # --------------------------------------------------------------------------------------------------------------
        # LAUNCH
        # --------------------------------------------------------------------------------------------------------------
        # p.solve_newton_2(mat, verbose=False, check=False)
        p.solve_newton_exact(mat, verbose=False, check=False)

        res_folder = "../../../res"

        def __plot(column: int):

            _, _, filenames = next(walk(res_folder))
            for time_step_index in range(len(time_steps)):
                for filename in filenames:
                    if "{}".format(time_step_index).zfill(6) in filename and "qdp" in filename:
                        hho_file_path = path.join(res_folder, filename)
                        with open(hho_file_path, "r") as hho_res_file:
                            fig, ax0d = plt.subplots(nrows=1, ncols=1)
                            c_hho = hho_res_file.readlines()
                            field_label = c_hho[0].split(",")[column]
                            number_of_points = len(c_hho)-1
                            eucli_d = displacement.euclidean_dimension
                            points = np.zeros((eucli_d, number_of_points), dtype=real)
                            field_vals = np.zeros((number_of_points, ), dtype=real)
                            for l_count, line in enumerate(c_hho[1:]):
                                x_coordinates = float(line.split(",")[0])
                                y_coordinates = float(line.split(",")[1])
                                field_value = float(line.split(",")[column])
                                points[0,l_count] += x_coordinates
                                points[1,l_count] += y_coordinates
                                field_vals[l_count] += field_value
                            x, y = points
                            colors = [(0, 0, 1), (0, 1, 1), (0, 1, 0), (1, 1, 0), (1, 0, 0)]
                            perso = LinearSegmentedColormap.from_list("perso", colors, N=1000)
                            vmin = min(field_vals[:])
                            vmax = max(field_vals[:])
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

        __plot(12)

        self.assertTrue(True)
