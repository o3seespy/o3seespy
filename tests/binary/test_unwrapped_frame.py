from collections import OrderedDict

import sfsimodels
import eqsig
import eqsig.duhamels

import numpy as np

from openseespy import opensees as op
from tests import extras as opc
from tests.conftest import TEST_DATA_DIR


def calc_yield_curvature(depth, eps_yield):
    """
    The yield curvature of a section from Priestley (Fig 4.15)

    :param depth:
    :param eps_yield:
    :return:
    """
    # TODO: get full validation of equation
    return 2.1 * eps_yield / depth


def elastic_bilin(ep1, ep2, eps_p2, en1=None, en2=None, eps_n2=None):
    if en1 is None:
        en1 = ep1
    if en2 is None:
        en2 = ep2
    if eps_n2 is None:
        eps_n2 = eps_p2

    return [ep1, ep2, eps_p2, en1, en2, eps_n2]


def get_inelastic_response(fb, asig, extra_time=0.0, xi=0.05, analysis_dt=0.001):
    """
    Run seismic analysis of a nonlinear FrameBuilding

    :param fb: FrameBuilding object
    :param asig: AccSignal object
    :param extra_time: float, additional analysis time after end of ground motion
    :param xi: damping ratio
    :param analysis_dt: time step to perform the analysis
    :return:
    """

    op.wipe()
    op.model('basic', '-ndm', 2, '-ndf', 3)  # 2 dimensions, 3 dof per node

    q_floor = 10000.  # kPa
    trib_width = fb.floor_length
    trib_mass_per_length = q_floor * trib_width / 9.8

    # Establish nodes and set mass based on trib area
    # Nodes named as: C<column-number>-S<storey-number>, first column starts at C1-S0 = ground level left
    nd = OrderedDict()
    col_xs = np.cumsum(fb.bay_lengths)
    col_xs = np.insert(col_xs, 0, 0)
    n_cols = len(col_xs)
    sto_ys = fb.heights
    sto_ys = np.insert(sto_ys, 0, 0)
    for cc in range(1, n_cols + 1):
        for ss in range(fb.n_storeys + 1):
            n_i = cc * 100 + ss
            nd["C%i-S%i" % (cc, ss)] = n_i
            op.node(n_i, col_xs[cc - 1], sto_ys[ss])

            if ss != 0:
                if cc == 1:
                    node_mass = trib_mass_per_length * fb.bay_lengths[0] / 2
                elif cc == n_cols:
                    node_mass = trib_mass_per_length * fb.bay_lengths[-1] / 2
                else:
                    node_mass = trib_mass_per_length * (fb.bay_lengths[cc - 2] + fb.bay_lengths[cc - 1] / 2)
                op.mass(n_i, node_mass)

    # Set all nodes on a storey to have the same displacement
    for ss in range(0, fb.n_storeys + 1):
        for cc in range(1, n_cols + 1):
            op.equalDOF(nd["C%i-S%i" % (1, ss)], nd["C%i-S%i" % (cc, ss)],  opc.X)

    # Fix all base nodes
    for cc in range(1, n_cols + 1):
        op.fix(nd["C%i-S%i" % (cc, 0)], opc.FIXED, opc.FIXED, opc.FIXED)

    # Coordinate transformation
    geo_tag = 1
    trans_args = []
    op.geomTransf("Linear", geo_tag, *[])

    l_hinge = fb.bay_lengths[0] * 0.1

    # Define material
    e_conc = 30.0e6
    i_beams = 0.4 * fb.beam_widths * fb.beam_depths ** 3 / 12
    i_columns = 0.5 * fb.column_widths * fb.column_depths ** 3 / 12
    a_beams = fb.beam_widths * fb.beam_depths
    a_columns = fb.column_widths * fb.column_depths
    ei_beams = e_conc * i_beams
    ei_columns = e_conc * i_columns
    eps_yield = 300.0e6 / 200e9
    phi_y_col = calc_yield_curvature(fb.column_depths, eps_yield)
    phi_y_beam = calc_yield_curvature(fb.beam_depths, eps_yield)

    # Define beams and columns

    md = OrderedDict()  # material dict
    sd = OrderedDict()  # section dict
    ed = OrderedDict()  # element dict
    # Columns named as: C<column-number>-S<storey-number>, first column starts at C1-S0 = ground floor left
    # Beams named as: B<bay-number>-S<storey-number>, first beam starts at B1-S1 = first storey left (foundation at S0)
    for ss in range(fb.n_storeys):

        # set columns
        for cc in range(1, fb.n_cols + 1):
            ele_i = cc * 100 + ss
            md["C%i-S%i" % (cc, ss)] = ele_i
            sd["C%i-S%i" % (cc, ss)] = ele_i
            ed["C%i-S%i" % (cc, ss)] = ele_i
            mat_props = elastic_bilin(ei_columns[ss][cc - 1], 0.05 * ei_columns[ss][cc - 1], phi_y_col[ss][cc - 1])
            #print(opc.ELASTIC_BILIN, ele_i, *mat_props)
            op.uniaxialMaterial(opc.ELASTIC_BILIN, ele_i, *mat_props)
            # op.uniaxialMaterial("Elastic", ele_i, ei_columns[ss][cc - 1])

            node_numbers = [nd["C%i-S%i" % (cc, ss)], nd["C%i-S%i" % (cc, ss + 1)]]
            op.element(opc.ELASTIC_BEAM_COLUMN, ele_i,
                       *node_numbers,
                       a_columns[ss - 1][cc - 1],
                       e_conc,
                       i_columns[ss - 1][cc - 1],
                       geo_tag
                       )

        # Set beams
        for bb in range(1, fb.n_bays + 1):
            ele_i = bb * 10000 + ss
            md["B%i-S%i" % (bb, ss)] = ele_i
            sd["B%i-S%i" % (bb, ss)] = ele_i
            ed["B%i-S%i" % (bb, ss)] = ele_i
            mat_props = elastic_bilin(ei_beams[ss][bb - 1], 0.05 * ei_beams[ss][bb - 1], phi_y_beam[ss][bb - 1])
            op.uniaxialMaterial(opc.ELASTIC_BILIN, ele_i, *mat_props)
            # op.uniaxialMaterial("Elastic", ele_i, ei_beams[ss][bb - 1])
            node_numbers = [nd["C%i-S%i" % (bb, ss + 1)], nd["C%i-S%i" % (bb + 1, ss + 1)]]
            print((opc.BEAM_WITH_HINGES, ele_i,
                       *node_numbers,
                        sd["B%i-S%i" % (bb, ss)], l_hinge,
                       sd["B%i-S%i" % (bb, ss)], l_hinge,
                       sd["B%i-S%i" % (bb, ss)], geo_tag
                       ))
            # Old definition
            # op.element(opc.BEAM_WITH_HINGES, ele_i,
            #  *[nd["C%i-S%i" % (bb, ss - 1)], nd["C%i-S%i" % (bb + 1, ss)]],
            #  sd["B%i-S%i" % (bb, ss)], l_hinge,
            #  sd["B%i-S%i" % (bb, ss)], l_hinge,
            #  e_conc,
            #  a_beams[ss - 1][bb - 1],
            #  i_beams[ss - 1][bb - 1], geo_tag
            #  )
            # New definition
            # op.element(opc.BEAM_WITH_HINGES, ele_i,
            #            *node_numbers,
            #             sd["B%i-S%i" % (bb, ss)], l_hinge,
            #            sd["B%i-S%i" % (bb, ss)], l_hinge,
            #            sd["B%i-S%i" % (bb, ss)], geo_tag  # TODO: make this elastic
            #
            #            )

            # Elastic definition
            op.element(opc.ELASTIC_BEAM_COLUMN, ele_i,
                       *node_numbers,
                       a_beams[ss - 1][bb - 1],
                       e_conc,
                       i_beams[ss - 1][bb - 1],
                       geo_tag
                       )

    # Define the dynamic analysis
    load_tag_dynamic = 1
    pattern_tag_dynamic = 1

    values = list(-1 * asig.values)  # should be negative
    op.timeSeries('Path', load_tag_dynamic, '-dt', asig.dt, '-values', *values)
    op.pattern('UniformExcitation', pattern_tag_dynamic, opc.X, '-accel', load_tag_dynamic)

    # set damping based on first eigen mode
    angular_freq = op.eigen('-fullGenLapack', 1)[0] ** 0.5
    if isinstance(angular_freq, complex):
        raise ValueError("Angular frequency is complex, issue with stiffness or mass")
    alpha_m = 0.0
    beta_k = 2 * xi / angular_freq
    beta_k_comm = 0.0
    beta_k_init = 0.0
    period = angular_freq / 2 / np.pi
    print("period: ", period)

    op.rayleigh(alpha_m, beta_k, beta_k_init, beta_k_comm)

    # Run the dynamic analysis

    op.wipeAnalysis()

    op.algorithm('Newton')
    op.system('SparseGeneral')
    op.numberer('RCM')
    op.constraints('Transformation')
    op.integrator('Newmark', 0.5, 0.25)
    op.analysis('Transient')
    #op.test("NormDispIncr", 1.0e-1, 2, 0)
    tol = 1.0e-10
    iter = 10
    op.test('EnergyIncr', tol, iter, 0, 2)  # TODO: make this test work
    analysis_time = (len(values) - 1) * asig.dt + extra_time
    outputs = {
        "time": [],
        "rel_disp": [],
        "rel_accel": [],
        "rel_vel": [],
        "force": []
    }
    print("Analysis starting")
    while op.getTime() < analysis_time:
        curr_time = op.getTime()
        op.analyze(1, analysis_dt)
        outputs["time"].append(curr_time)
        outputs["rel_disp"].append(op.nodeDisp(nd["C%i-S%i" % (1, fb.n_storeys)], opc.X))
        outputs["rel_vel"].append(op.nodeVel(nd["C%i-S%i" % (1, fb.n_storeys)], opc.X))
        outputs["rel_accel"].append(op.nodeAccel(nd["C%i-S%i" % (1, fb.n_storeys)], opc.X))
        op.reactions()
        react = 0
        for cc in range(1, fb.n_cols):
            react += -op.nodeReaction(nd["C%i-S%i" % (cc, 0)], opc.X)
        outputs["force"].append(react)  # Should be negative since diff node
    op.wipe()
    for item in outputs:
        outputs[item] = np.array(outputs[item])

    return outputs

#
# def load_frame_building_sample_data():
#     """
#     Sample data for the FrameBuilding object
#
#     :param fb:
#     :return:
#     """
#     number_of_storeys = 6
#     interstorey_height = 3.4  # m
#     masses = 40.0e3  # kg
#     n_bays = 3
#
#     fb = sfsimodels.FrameBuilding2D(number_of_storeys, n_bays)
#     fb.interstorey_heights = interstorey_height * np.ones(number_of_storeys)
#     fb.floor_length = 18.0  # m
#     fb.floor_width = 16.0  # m
#     fb.storey_masses = masses * np.ones(number_of_storeys)  # kg
#
#     fb.bay_lengths = [6., 6.0, 6.0]
#     fb.set_beam_prop("depth", [0.5, 0.5, 0.5], repeat="up")
#     fb.set_beam_prop("width", [0.4, 0.4, 0.4], repeat="up")
#     fb.set_column_prop("width", [0.5, 0.5, 0.5, 0.5], repeat="up")
#     fb.set_column_prop("depth", [0.5, 0.5, 0.5, 0.5], repeat="up")
#     fb.n_seismic_frames = 3
#     fb.n_gravity_frames = 0
#     return fb


def load_small_frame_building_sample_data():
    """
    Sample data for the FrameBuilding object

    :param fb:
    :return:
    """
    number_of_storeys = 1
    interstorey_height = 3.4  # m
    masses = 40.0e3  # kg
    n_bays = 1

    fb = sfsimodels.FrameBuilding2D(number_of_storeys, n_bays)
    fb.interstorey_heights = interstorey_height * np.ones(number_of_storeys)
    fb.floor_length = 18.0  # m
    fb.floor_width = 16.0  # m
    fb.storey_masses = masses * np.ones(number_of_storeys)  # kg

    fb.bay_lengths = [6.]
    fb.set_beam_prop("depth", [0.5], repeat="up")
    fb.set_beam_prop("width", [0.4], repeat="up")
    fb.set_column_prop("width", [0.5, 0.5], repeat="up")
    fb.set_column_prop("depth", [0.5, 0.5], repeat="up")
    return fb


def plot_response():
    import matplotlib.pyplot as plt

    record_path = TEST_DATA_DIR
    record_filename = 'test_motion_dt0p01.txt'
    dt = 0.01
    rec = np.loadtxt(record_path + record_filename)
    acc_signal = eqsig.AccSignal(rec, dt)
    xi = 0.05
    time = acc_signal.time

    # frame = load_frame_building_sample_data()
    frame = load_small_frame_building_sample_data()
    print("Building loaded")

    outputs = get_inelastic_response(frame, acc_signal, xi=xi, extra_time=0)
    print("Analysis complete")
    acc_opensees = np.interp(time, outputs["time"], outputs["rel_accel"]) - rec
    ux_opensees = np.interp(time, outputs["time"], outputs["rel_disp"])
    plt.plot(time, ux_opensees, label="BB")

    periods = np.array([2 * 0.7147])
    resp_u, resp_v, resp_a = eqsig.duhamels.response_series(motion=rec, dt=dt, periods=periods, xi=xi)
    plt.plot(acc_signal.time, resp_u[0], label="Duhamels")
    print(np.sum(np.abs(resp_u[0])), )
    assert np.isclose(np.sum(np.abs(ux_opensees)), 113.09023), np.sum(np.abs(ux_opensees))
    plt.legend()
    plt.show()
    print("Complete")


def test_small_frame_dynamic():
    record_path = TEST_DATA_DIR
    record_filename = 'test_motion_dt0p01.txt'
    dt = 0.01
    rec = np.loadtxt(record_path + record_filename)
    acc_signal = eqsig.AccSignal(rec, dt)
    xi = 0.05
    time = acc_signal.time

    frame = load_small_frame_building_sample_data()

    outputs = get_inelastic_response(frame, acc_signal, xi=xi, extra_time=0)
    ux_opensees = np.interp(time, outputs["time"], outputs["rel_disp"])
    assert np.isclose(np.sum(np.abs(ux_opensees)), 113.09023), np.sum(np.abs(ux_opensees))


if __name__ == '__main__':

    plot_response()
