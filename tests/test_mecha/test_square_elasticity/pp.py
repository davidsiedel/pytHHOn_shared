from unittest import TestCase

from numpy import ndarray

from pythhon.pbbb.problem import Problem
from pythhon.pbbb.boundary_condition import BoundaryCondition
from pythhon.pbbb.load import Load
from pythhon.pbbb.field import Field
from pythhon.fem.element.finite_element import FiniteElement
# from pythhon.parameters import *
# from pp.post_processing import *

from mgis import behaviour as mgis_bv

from pythhon.pbbb.material import Material

import matplotlib.pyplot as plt
from os import walk, path
import numpy as np

def plot_data_2(
        mtest_file_path: str,
        hho_res_dir_path: str,
        number_of_time_steps: int,
        m_x_inedx: int,
        m_y_index: int,
        d_x_inedx: int,
        d_y_inedx: int,
):
    coef = 1.0e6
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
    # plt.ylim(-0.01* 1.e-8, 0.009 * 1.e-7)
    # plt.ylim(-50.0, 1000.0)  # SHEAR
    # plt.xlim(-0.05, 0.05)
    plt.grid()
    plt.show()
    return
