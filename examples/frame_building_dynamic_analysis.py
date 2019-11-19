from collections import OrderedDict

import matplotlib.pyplot as plt
import sfsimodels
import eqsig

import numpy as np

from openseespy import opensees as opy
import o3seespy as o3


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

    Parameters
    ----------
    fb: sfsimodels.Frame2DBuilding object
    asig: eqsig.AccSignal object
    extra_time
    xi
    analysis_dt

    Returns
    -------

    """
    osi = o3.OpenseesInstance(dimensions=2)

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
            nd["C{0}-S{1}".format(cc, ss)] = o3.node.Node(osi, col_xs[cc - 1], sto_ys[ss])

            if ss != 0:
                if cc == 1:
                    node_mass = trib_mass_per_length * fb.bay_lengths[0] / 2
                elif cc == n_cols:
                    node_mass = trib_mass_per_length * fb.bay_lengths[-1] / 2
                else:
                    node_mass = trib_mass_per_length * (fb.bay_lengths[cc - 2] + fb.bay_lengths[cc - 1] / 2)
                o3.set_node_mass(nd["C{0}-S{1}".format(cc, ss)], node_mass, 0., 0.)

    # Set all nodes on a storey to have the same displacement
    for ss in range(0, fb.n_storeys + 1):
        for cc in range(1, n_cols + 1):
            o3.set_equal_dof(nd["C{0}-S{1}".format(1, ss)], nd["C{0}-S{1}".format(cc, ss)], o3.cc.X)

    # Fix all base nodes
    for cc in range(1, n_cols + 1):
        opy.fix(nd["C%i-S%i" % (cc, 0)].tag, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED)

    # Coordinate transformation
    transf = o3.transformation.Linear(osi, [])

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
    phi_y_beam = calc_yield_curvature(fb.beam_depths, eps_yield) * 10  # TODO: re-evaluate

    # Define beams and columns
    # Columns named as: C<column-number>-S<storey-number>, first column starts at C1-S0 = ground floor left
    # Beams named as: B<bay-number>-S<storey-number>, first beam starts at B1-S1 = first storey left (foundation at S0)

    md = OrderedDict()  # material dict
    sd = OrderedDict()  # section dict
    ed = OrderedDict()  # element dict

    for ss in range(fb.n_storeys):

        # set columns
        for cc in range(1, fb.n_cols + 1):
            lp_i = 0.4
            lp_j = 0.4  # plastic hinge length
            ele_str = "C{0}-S{1}S{2}".format(cc, ss, ss + 1)

            top_sect = o3.section.Elastic2D(osi, e_conc, a_columns[ss][cc - 1], i_columns[ss][cc - 1])
            bot_sect = o3.section.Elastic2D(osi, e_conc, a_columns[ss][cc - 1], i_columns[ss][cc - 1])
            centre_sect = o3.section.Elastic2D(osi, e_conc, a_columns[ss][cc - 1], i_columns[ss][cc - 1])
            sd[ele_str + "T"] = top_sect
            sd[ele_str + "B"] = bot_sect
            sd[ele_str + "C"] = centre_sect

            integ = o3.beam_integration.HingeMidpoint(osi, bot_sect, lp_i, top_sect, lp_j, centre_sect)

            left_node = nd["C%i-S%i" % (cc, ss)]
            right_node = nd["C%i-S%i" % (cc, ss + 1)]
            ed[ele_str] = o3.element.ForceBeamColumn(osi, left_node, right_node, transf, integ)

        # Set beams
        for bb in range(1, fb.n_bays + 1):
            lp_i = 0.5
            lp_j = 0.5
            ele_str = "C{0}C{1}-S{2}".format(bb - 1, bb, ss)

            mat = o3.uniaxial_material.ElasticBiLinear(osi, ei_beams[ss][bb - 1], 0.05 * ei_beams[ss][bb - 1], phi_y_beam[ss][bb - 1])
            md[ele_str] = mat
            left_sect = o3.section.Uniaxial(osi, mat, quantity=o3.cc.M_Z)
            right_sect = o3.section.Uniaxial(osi, mat, quantity=o3.cc.M_Z)
            centre_sect = o3.section.Elastic2D(osi, e_conc, a_beams[ss][bb - 1], i_beams[ss][bb - 1])
            integ = o3.beam_integration.HingeMidpoint(osi, left_sect, lp_i, right_sect, lp_j, centre_sect)

            left_node = nd["C%i-S%i" % (bb, ss + 1)]
            right_node = nd["C%i-S%i" % (bb + 1, ss + 1)]
            ed[ele_str] = o3.element.ForceBeamColumn(osi, left_node, right_node, transf, integ)

    # Define the dynamic analysis
    load_tag_dynamic = 1
    pattern_tag_dynamic = 1

    values = list(-1 * asig.values)  # should be negative
    opy.timeSeries('Path', load_tag_dynamic, '-dt', asig.dt, '-values', *values)
    opy.pattern('UniformExcitation', pattern_tag_dynamic, o3.cc.X, '-accel', load_tag_dynamic)

    # set damping based on first eigen mode
    angular_freq = opy.eigen('-fullGenLapack', 1) ** 0.5
    if isinstance(angular_freq, complex):
        raise ValueError("Angular frequency is complex, issue with stiffness or mass")
    beta_k = 2 * xi / angular_freq
    o3.rayleigh.Rayleigh(osi, alpha_m=0.0, beta_k=beta_k, beta_k_init=0.0, beta_k_comm=0.0)

    # Run the dynamic analysis

    opy.wipeAnalysis()

    o3.algorithm.Newton(osi)
    opy.system('SparseGeneral')
    opy.numberer('RCM')
    opy.constraints('Transformation')
    opy.integrator('Newmark', 0.5, 0.25)
    opy.analysis('Transient')
    #op.test("NormDispIncr", 1.0e-1, 2, 0)
    tol = 1.0e-4
    iter = 4
    opy.test('EnergyIncr', tol, iter, 0, 2)
    analysis_time = (len(values) - 1) * asig.dt + extra_time
    outputs = {
        "time": [],
        "rel_disp": [],
        "rel_accel": [],
        "rel_vel": [],
        "force": [],
        "ele_mom": [],
        "ele_curve": [],
    }
    print("Analysis starting")
    opy.recorder('Element', '-file', 'ele_out.txt', '-time', '-ele', 1, 'force')
    while opy.getTime() < analysis_time:
        curr_time = opy.getTime()
        opy.analyze(1, analysis_dt)
        outputs["time"].append(curr_time)
        outputs["rel_disp"].append(opy.nodeDisp(nd["C%i-S%i" % (1, fb.n_storeys)].tag, o3.cc.X))
        outputs["rel_vel"].append(opy.nodeVel(nd["C%i-S%i" % (1, fb.n_storeys)].tag, o3.cc.X))
        outputs["rel_accel"].append(opy.nodeAccel(nd["C%i-S%i" % (1, fb.n_storeys)].tag, o3.cc.X))
        # outputs['ele_mom'].append(opy.eleResponse('-ele', [ed['B%i-S%i' % (1, 0)], 'basicForce']))
        opy.reactions()
        react = 0
        for cc in range(1, fb.n_cols):
            react += -opy.nodeReaction(nd["C%i-S%i" % (cc, 0)].tag, o3.cc.X)
        outputs["force"].append(react)  # Should be negative since diff node
    opy.wipe()
    for item in outputs:
        outputs[item] = np.array(outputs[item])

    return outputs


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


if __name__ == '__main__':
    from tests import conftest
    record_filename = 'test_motion_dt0p01.txt'
    motion_step = 0.01
    rec = np.loadtxt(conftest.TEST_DATA_DIR + record_filename)

    xi = 0.05

    acc_signal = eqsig.AccSignal(rec, motion_step)
    time = acc_signal.time

    # frame = load_frame_building_sample_data()
    frame = load_small_frame_building_sample_data()
    print("Building loaded")

    outputs = get_inelastic_response(frame, acc_signal, xi=xi, extra_time=0)
    print("Analysis complete")
    acc_opensees = np.interp(time, outputs["time"], outputs["rel_accel"]) - rec
    ux_opensees = np.interp(time, outputs["time"], outputs["rel_disp"])
    plt.plot(ux_opensees)
    plt.show()
    print("Complete")
