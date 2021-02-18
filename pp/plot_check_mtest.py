import matplotlib.pyplot as plt


def plot_data(
    mtest_file_path: str,
    hho_file_path: str,
    m_x_inedx: int,
    m_y_index: int,
    d_x_min_inedx: int,
    d_x_max_inedx: int,
    d_y_min_index: int,
    d_y_max_index: int,
):
    with open(mtest_file_path, "r") as mres_file:
        c = mres_file.readlines()
        labels = []
        m_x_data = []
        m_y_data = []
        for line in c:
            if line[0] == "#":
                alpha = line.split("column:")
                labels.append(alpha[1].replace("\n", ""))
        coef = 10.0e6
        line_start = len(labels)
        for time, line in enumerate(c[line_start:]):
            m_x_data.append(float(line.split(" ")[m_x_inedx]))
            m_y_data.append(float(line.split(" ")[m_y_index]) / coef)
        m_x_label = labels[m_x_inedx]
        m_y_label = labels[m_y_index]
    # # with open("/home/dsiedel/Projects/pythhon/behaviour/testfront/small_strain_isotropic_linear_kinematic_hardening.csv", "r") as mres_file:
    # with open("/home/dsiedel/Projects/pythhon/behaviour/testfront/finite_strain_isotropic_linear_hardening.csv", "r") as mres_file:
    #     c = mres_file.readlines()
    #     labels = c[0].split(" ")
    #     m_x_data = []
    #     m_y_data = []
    #     # 0 -> eps_xx
    #     # 1 -> eps_yy
    #     # 2 -> eps_zz
    #     # 3 -> eps_xy
    #     # 4 -> sig_xx
    #     # 5 -> sig_yy
    #     # 6 -> sig_zz
    #     # 7 -> sig_xy
    #     for time, line in enumerate(c[1:]):
    #         m_x_data.append(float(line.split(" ")[m_x_inedx]))
    #         m_y_data.append(float(line.split(" ")[m_y_index])/10.e6)
    #         m_x_label = labels[m_x_inedx]
    #         m_y_label = labels[m_y_index]
    with open(hho_file_path, "r") as d_res_file:
        c = d_res_file.readlines()
        labels = c[0].split(" ")
        d_x_min_data = []
        d_y_min_data = []
        d_x_max_data = []
        d_y_max_data = []
        # 0 -> eps_xx
        # 1 -> eps_yy
        # 2 -> eps_zz
        # 3 -> eps_xy
        # 4 -> sig_xx
        # 5 -> sig_yy
        # 6 -> sig_zz
        # 7 -> sig_xy
        for time, line in enumerate(c[1:]):
            d_x_min_data.append(float(line.split(" ")[d_x_min_inedx]))
            d_y_min_data.append(float(line.split(" ")[d_y_min_index]) / 10.0e6)
            d_x_max_data.append(float(line.split(" ")[d_x_max_inedx]))
            d_y_max_data.append(float(line.split(" ")[d_y_max_index]) / 10.0e6)
    plt.plot(m_x_data, m_y_data, color="blue", label="MTEST")
    plt.plot(d_x_min_data, d_y_min_data, color="green", linestyle="--", label="HHO_MIN")
    plt.plot(d_x_max_data, d_y_max_data, color="purple", linestyle="--", label="HHO_MAX")
    plt.xlabel(m_x_label)
    plt.ylabel(m_y_label + " [MPA]")
    plt.title("MTEST HHO COMPARISON CYCLIC LOADING ISOTROPIC LINEAR HARDENING")
    plt.legend()
    # plt.ylim(-0.35, 0.35) #TRACTION
    plt.ylim(-200.0, 200.0)  # SHEAR
    # plt.xlim(-0.05, 0.05)
    plt.grid()
    plt.show()
    # eps_xx = line.split(" ")[0]
    # eps_yy = line.split(" ")[1]
    # eps_zz = line.split(" ")[2]
    # eps_xy = line.split(" ")[3]
    # sig_xx = line.split(" ")[4]
    # sig_yy = line.split(" ")[5]
    # sig_zz = line.split(" ")[6]
    # sig_xy = line.split(" ")[7]

# mtest_file_path = "/home/dsiedel/Projects/pythhon/behaviour/testfront/elasticity.res"
# hho_file_path = "/home/dsiedel/Projects/pythhon/res/res.csv"
# # TRACTION COMPRESSION
# plot_data(mtest_file_path, hho_file_path, 1, 5, 1, 2, 9, 10) # EPS_XX - SIG_XX
# plot_data(mtest_file_path, hho_file_path, 1, 6, 1, 2, 11, 12) # EPS_XX - SIG_YY
# plot_data(mtest_file_path, hho_file_path, 1, 7, 1, 2, 13, 14) # EPS_XX - SIG_ZZ
# plot_data(mtest_file_path, hho_file_path, 1, 8, 1, 2, 15, 16) # EPS_XX - SIG_XY

# SHEAR
# plot_data(4, 8, 7, 8, 15, 16) # EPS_XY - SIG_XY
# plot_data(4, 5, 7, 8, 9, 10) # EPS_XX - SIG_YY
# plot_data(4, 6, 7, 8, 11, 12) # EPS_XX - SIG_ZZ
# plot_data(4, 7, 7, 8, 13, 14) # EPS_XX - SIG_XY

# # TRACTION COMPRESSION FINITE STRAIN
# plot_data(1, 6, 1, 2, 9, 10)  # EPS_XX - SIG_XX
# plot_data(1, 7, 1, 2, 11, 12) # EPS_XX - SIG_YY
# plot_data(1, 8, 1, 2, 13, 14) # EPS_XX - SIG_ZZ
# plot_data(1, 9, 1, 2, 15, 16) # EPS_XX - SIG_XY
