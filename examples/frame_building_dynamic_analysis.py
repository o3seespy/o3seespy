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
    osi = o3.OpenSeesInstance(ndm=2)

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
            nd[f"C{cc}-S{ss}"] = o3.node.Node(osi, col_xs[cc - 1], sto_ys[ss])

            if ss != 0:
                if cc == 1:
                    node_mass = trib_mass_per_length * fb.bay_lengths[0] / 2
                elif cc == n_cols:
                    node_mass = trib_mass_per_length * fb.bay_lengths[-1] / 2
                else:
                    node_mass = trib_mass_per_length * (fb.bay_lengths[cc - 2] + fb.bay_lengths[cc - 1] / 2)
                o3.set_node_mass(osi, nd[f"C{cc}-S{ss}"], node_mass, 0., 0.)

    # Set all nodes on a storey to have the same displacement
    for ss in range(0, fb.n_storeys + 1):
        for cc in range(1, n_cols + 1):
            o3.set_equal_dof(osi, nd[f"C1-S{ss}"], nd[f"C{cc}-S{ss}"], o3.cc.X)

    # Fix all base nodes
    for cc in range(1, n_cols + 1):
        o3.Fix3DOF(osi, nd[f"C{cc}-S0"], o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED)

    # Coordinate transformation
    transf = o3.geom_transf.Linear2D(osi, [])

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
            ele_str = f"C{cc}-S{ss}S{ss+1}"

            top_sect = o3.section.Elastic2D(osi, e_conc, a_columns[ss][cc - 1], i_columns[ss][cc - 1])
            bot_sect = o3.section.Elastic2D(osi, e_conc, a_columns[ss][cc - 1], i_columns[ss][cc - 1])
            centre_sect = o3.section.Elastic2D(osi, e_conc, a_columns[ss][cc - 1], i_columns[ss][cc - 1])
            sd[ele_str + "T"] = top_sect
            sd[ele_str + "B"] = bot_sect
            sd[ele_str + "C"] = centre_sect

            integ = o3.beam_integration.HingeMidpoint(osi, bot_sect, lp_i, top_sect, lp_j, centre_sect)

            bot_node = nd[f"C{cc}-S{ss}"]
            top_node = nd[f"C{cc}-S{ss+1}"]
            ed[ele_str] = o3.element.ForceBeamColumn(osi, [bot_node, top_node], transf, integ)

        # Set beams
        for bb in range(1, fb.n_bays + 1):
            lp_i = 0.5
            lp_j = 0.5
            ele_str = f"C{bb-1}C{bb}-S{ss}"

            mat = o3.uniaxial_material.ElasticBilin(osi, ei_beams[ss][bb - 1], 0.05 * ei_beams[ss][bb - 1], phi_y_beam[ss][bb - 1])
            md[ele_str] = mat
            left_sect = o3.section.Uniaxial(osi, mat, quantity=o3.cc.M_Z)
            right_sect = o3.section.Uniaxial(osi, mat, quantity=o3.cc.M_Z)
            centre_sect = o3.section.Elastic2D(osi, e_conc, a_beams[ss][bb - 1], i_beams[ss][bb - 1])
            integ = o3.beam_integration.HingeMidpoint(osi, left_sect, lp_i, right_sect, lp_j, centre_sect)

            left_node = nd[f"C{bb}-S{ss+1}"]
            right_node = nd[f"C{bb+1}-S{ss+1}"]
            ed[ele_str] = o3.element.ForceBeamColumn(osi, [left_node, right_node], transf, integ)

    # Define the dynamic analysis
    a_series = o3.time_series.Path(osi, dt=asig.dt, values=-1 * asig.values)  # should be negative
    o3.pattern.UniformExcitation(osi, dir=o3.cc.X, accel_series=a_series)

    # set damping based on first eigen mode
    angular_freq = o3.get_eigen(osi, solver='fullGenLapack', n=1) ** 0.5
    if isinstance(angular_freq, complex):
        raise ValueError("Angular frequency is complex, issue with stiffness or mass")
    beta_k = 2 * xi / angular_freq
    o3.rayleigh.Rayleigh(osi, alpha_m=0.0, beta_k=beta_k, beta_k_init=0.0, beta_k_comm=0.0)

    # Run the dynamic analysis

    o3.wipe_analysis(osi)

    o3.algorithm.Newton(osi)
    o3.system.SparseGeneral(osi)
    o3.numberer.RCM(osi)
    o3.constraints.Transformation(osi)
    o3.integrator.Newmark(osi, 0.5, 0.25)
    o3.analysis.Transient(osi)

    tol = 1.0e-4
    iter = 4
    o3.test_check.EnergyIncr(osi, tol, iter, 0, 2)
    analysis_time = (len(asig.values) - 1) * asig.dt + extra_time
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
    while o3.get_time(osi) < analysis_time:
        curr_time = opy.getTime()
        o3.analyze(osi, 1, analysis_dt)
        outputs["time"].append(curr_time)
        outputs["rel_disp"].append(o3.get_node_disp(osi, nd["C%i-S%i" % (1, fb.n_storeys)], o3.cc.X))
        outputs["rel_vel"].append(o3.get_node_vel(osi, nd["C%i-S%i" % (1, fb.n_storeys)], o3.cc.X))
        outputs["rel_accel"].append(o3.get_node_accel(osi, nd["C%i-S%i" % (1, fb.n_storeys)], o3.cc.X))
        # outputs['ele_mom'].append(opy.eleResponse('-ele', [ed['B%i-S%i' % (1, 0)], 'basicForce']))
        o3.gen_reactions(osi)
        react = 0
        for cc in range(1, fb.n_cols):
            react += -o3.get_node_reaction(osi, nd["C%i-S%i" % (cc, 0)], o3.cc.X)
        outputs["force"].append(react)  # Should be negative since diff node
    o3.wipe(osi)
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
