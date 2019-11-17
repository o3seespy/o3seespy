import eqsig
from collections import OrderedDict
import numpy as np
import sfsimodels as sm
import openseespy.opensees as opy
import o3seespy as o3

# for linear analysis comparison
import liquepy as lq
import pysra
from bwplot import cbox


def site_response(sp, asig):
    """
    Run seismic analysis of a soil profile - example based on:
    http://opensees.berkeley.edu/wiki/index.php/Site_Response_Analysis_of_a_Layered_Soil_Column_(Total_Stress_Analysis)

    :param sp: sfsimodels.SoilProfile object
        A soil profile
    :param asig: eqsig.AccSignal object
        An acceleration signal
    :return:
    """
    # TODO: THIS IS NOT WORKING CORRECTLY, RESPONSE IS WRONG
    osi = o3.OpenseesInstance(dimensions=2, node_dofs=2, state=3)
    assert isinstance(sp, sm.SoilProfile)
    sp.gen_split(props=['shear_vel', 'unit_mass', 'cohesion', 'phi', 'bulk_mod', 'poissons_ratio', 'strain_peak'])
    thicknesses = sp.split["thickness"]
    n_node_rows = len(thicknesses) + 1
    node_depths = np.cumsum(sp.split["thickness"])
    node_depths = np.insert(node_depths, 0, 0)
    ele_depths = (node_depths[1:] + node_depths[:-1]) / 2
    shear_vels = sp.split["shear_vel"]
    unit_masses = sp.split["unit_mass"] / 1e3
    g_mods = unit_masses * shear_vels ** 2 / 1e3
    poissons_ratio = sp.split['poissons_ratio']
    youngs_mods = 2 * g_mods * (1 - poissons_ratio)
    bulk_mods = youngs_mods / (3 * (1 - 2 * poissons_ratio))
    bulk_mods = sp.split['bulk_mod'] / 1e3

    ref_pressure = 80.0
    cohesions = sp.split['cohesion']
    phis = sp.split['phi']
    strain_peaks = sp.split['strain_peak']
    grav = 9.81
    damping = 0.02
    omega_1 = 2 * np.pi * 0.2
    omega_2 = 2 * np.pi * 20
    a0 = 2 * damping * omega_1 * omega_2 / (omega_1 + omega_2)
    a1 = 2 * damping / (omega_1 + omega_2)

    newmark_gamma = 0.5
    newmark_beta = 0.25

    ele_width = min(thicknesses)
    total_soil_nodes = len(thicknesses) * 2 + 2

    # Define nodes and set boundary conditions for simple shear deformation
    # Start at top and build down?
    nd = OrderedDict()
    nd["R0L"] = o3.node.Node(osi, 0, 0)  # row 0 left
    nd["R0R"] = o3.node.Node(osi, ele_width, 0)
    for i in range(1, n_node_rows):
        # Establish left and right nodes
        nd["R{0}L".format(i)] = o3.node.Node(osi, 0, -node_depths[i])
        nd["R{0}R".format(i)] = o3.node.Node(osi, ele_width, -node_depths[i], x_mass=0.0001) # TODO: why is mass needed for stability?
        # set x and y dofs equal for left and right nodes
        if 1 == 0:
            if i != n_node_rows - 1:  # TODO: why not
                o3.EqualDOF(osi, nd["R{0}L".format(i)], nd["R{0}R".format(i)], [o3.cc.X, o3.cc.Y])

    # Fix base nodes
    o3.Fix(osi, nd["R{0}L".format(n_node_rows - 1)], o3.cc.FREE, o3.cc.FIXED, o3.cc.FREE)
    o3.Fix(osi, nd["R{0}R".format(n_node_rows - 1)], o3.cc.FREE, o3.cc.FIXED, o3.cc.FREE)

    # Define dashpot nodes
    dashpot_node_l = o3.node.Node(osi, 0, -node_depths[-1])
    dashpot_node_2 = o3.node.Node(osi, 0, -node_depths[-1])
    o3.Fix(osi, dashpot_node_l,  o3.cc.FIXED, o3.cc.FIXED, o3.cc.FREE)
    o3.Fix(osi, dashpot_node_2, o3.cc.FREE, o3.cc.FIXED, o3.cc.FREE)

    # define equal DOF for dashpot and soil base nodes
    o3.EqualDOF(osi, nd["R{0}L".format(n_node_rows - 1)], nd["R{0}R".format(n_node_rows - 1)], [o3.cc.X])
    o3.EqualDOF(osi, nd["R{0}L".format(n_node_rows - 1)], dashpot_node_2, [o3.cc.X])

    # define materials
    ele_thick = 1.0  # m
    soil_mats = []
    strains = np.logspace(-6, -0.5, 16)
    ref_strain = 0.005
    rats = 1. / (1 + (strains / ref_strain) ** 0.91)
    for i in range(len(thicknesses)):
        # mat = o3.nd_material.PressureIndependMultiYield(osi, 2, unit_masses[i], g_mods[i],
        #                                                  bulk_mods[i], cohesions[i], strain_peaks[i],
        #                                                  phis[i], press_depend_coe=0.0, no_yield_surf=16,
        #                                                  strains=strains, ratios=rats)
        mat = o3.nd_material.ElasticIsotropic(osi, youngs_mods[i], poissons_ratio[i])
        soil_mats.append(mat)

        # def element
        nodes = [

            nd["R{0}L".format(i + 1)],
            nd["R{0}R".format(i + 1)],
            nd["R{0}R".format(i)],
            nd["R{0}L".format(i)]
        ]
        ele = o3.element.Quad(osi, nodes, ele_thick, o3.cc.PLANE_STRAIN, mat, b2=grav * unit_masses[i])

    # define material and element for viscous dampers
    dashpot_c = ele_width * unit_masses[-1] * shear_vels[-1]
    dashpot_mat = o3.uniaxial_material.Viscous(osi, dashpot_c, alpha=1.)
    dashpot_ele = o3.element.ZeroLength(osi, dashpot_node_l, dashpot_node_2, mat_x=dashpot_mat)

    o3.constraint.Transformation(osi)
    o3.test_check.NormDispIncr(osi, tol=1.0e-5, max_iter=30, p_flag=0)
    o3.algorithm.Newton(osi)
    o3.numberer.RCM(osi)
    o3.system.ProfileSPD(osi)
    o3.integrator.Newmark(osi, newmark_gamma, newmark_beta)

    o3.analysis.Transient(osi)
    # opy.analyze(10, 5.0e3)
    o3.analyze(osi, 10, 500.)

    for i in range(len(soil_mats)):
        o3.update_material_stage(osi, soil_mats[i], 1)
        # opy.updateMaterialStage('-material', soil_mats[i].tag, '-stage', 1)
    o3.analyze(osi, 10, 500.)
    # reset time and analysis
    o3.set_time(osi, 0.0)
    opy.setTime(0.0)
    o3.wipe_analysis(osi)

    # Define the dynamic analysis
    load_tag_dynamic = 1
    pattern_tag_dynamic = 1

    values = list(-1 * asig.values)  # should be negative
    opy.timeSeries('Path', load_tag_dynamic, '-dt', asig.dt, '-values', *values)
    opy.pattern('UniformExcitation', pattern_tag_dynamic, o3.cc.X, '-accel', load_tag_dynamic)

    # set damping based on first eigen mode
    xi = 0.01
    # angular_freq = opy.eigen('-fullGenLapack', 1) ** 0.5
    angular_freq = 0.5
    beta_k = 2 * xi / angular_freq
    o3.rayleigh.Rayleigh(osi, alpha_m=0.0, beta_k=beta_k, beta_k_init=0.0, beta_k_comm=0.0)

    # Run the dynamic analysis
    o3.algorithm.Newton(osi)
    o3.system.SparseGeneral(osi)
    o3.numberer.RCM(osi)
    o3.constraint.Transformation(osi)
    o3.integrator.Newmark(osi, newmark_gamma, newmark_beta)
    o3.analysis.Transient(osi)

    o3.test_check.EnergyIncr(osi, tol=1.0e-10, max_iter=10)
    analysis_time = (len(values) - 1) * asig.dt
    analysis_dt = 0.01

    o3.recorder.NodeToFile(osi, 'sample_out.txt', node=nd["R0L"], dofs=[o3.cc.X], res_type='accel')
    na = o3.recorder.NodeToArrayCache(osi, node=nd["R0L"], dofs=[o3.cc.X], res_type='accel')

    while opy.getTime() < analysis_time:
        print(opy.getTime())
        opy.analyze(1, analysis_dt)
    opy.wipe()
    outputs = {
        "time": np.arange(0, analysis_time, analysis_dt),
        "rel_disp": [],
        "rel_accel": na.collect(),
        "rel_vel": [],
        "force": []
    }

    return outputs


def run_pysra(soil_profile, asig, odepths):
    pysra_profile = lq.sra.sm_profile_to_pysra(soil_profile, d_inc=[0.5] * soil_profile.n_layers)
    pysra_m = pysra.motion.TimeSeriesMotion(filename=asig.label, description=None, time_step=asig.dt,
                                      accels=-asig.values / 9.8)  # Should be input as g

    calc = pysra.propagation.LinearElasticCalculator()

    od = {}
    outs = []
    for i, depth in enumerate(odepths):
        od["ACCX_d%i" % i] = len(outs)
        outs.append(pysra.output.AccelerationTSOutput(pysra.output.OutputLocation('within', depth=depth)))
    outputs = pysra.output.OutputCollection(outs)
    calc(pysra_m, pysra_profile, pysra_profile.location('outcrop', depth=soil_profile.height))
    outputs(calc)

    out_series = {}
    for item in od:
        out_series[item] = outputs[od[item]].values[:asig.npts] * 9.8
    return out_series


def run():
    sl = sm.Soil()
    # vs = 250.
    vs = 100.
    unit_mass = 1700.0
    sl.g_mod = vs ** 2 * unit_mass
    sl.poissons_ratio = 0.0
    sl.cohesion = 95.0e3
    sl.phi = 0.0
    sl.unit_dry_weight = unit_mass * 9.8
    sl.strain_peak = 0.1  # set additional parameter required for PIMY model
    sl.xi = 0.03  # for linear analysis
    assert np.isclose(vs, sl.get_shear_vel(saturated=False))
    soil_profile = sm.SoilProfile()
    soil_profile.add_layer(0, sl)
    sl = sm.Soil()
    # vs = 250.
    vs = 400.
    unit_mass = 1700.0
    sl.g_mod = vs ** 2 * unit_mass
    sl.poissons_ratio = 0.0
    sl.cohesion = 395.0e3
    sl.phi = 0.0
    sl.unit_dry_weight = unit_mass * 9.8
    sl.strain_peak = 0.1  # set additional parameter required for PIMY model
    sl.xi = 0.03  # for linear analysis
    soil_profile.add_layer(6, sl)
    soil_profile.height = 12.0
    from tests.conftest import TEST_DATA_DIR

    record_path = TEST_DATA_DIR
    record_filename = 'test_motion_dt0p01.txt'
    dt = 0.01
    rec = np.loadtxt(record_path + record_filename) / 10
    acc_signal = eqsig.AccSignal(rec, dt)

    import matplotlib.pyplot as plt
    bf, sps = plt.subplots(nrows=3)

    # linear analysis with pysra
    od = run_pysra(soil_profile, acc_signal, odepths=np.array([0.0, 2.0]))
    pysra_sig = eqsig.AccSignal(od['ACCX_d0'], acc_signal.dt)

    outputs = site_response(soil_profile, acc_signal)
    resp_dt = outputs['time'][2] - outputs['time'][1]
    surf_sig = eqsig.AccSignal(outputs['rel_accel'], resp_dt)  # TODO: need to get absolute acceleration

    sps[0].plot(acc_signal.time, acc_signal.values, c='k')
    sps[0].plot(surf_sig.time, surf_sig.values, c=cbox(0))
    sps[0].plot(acc_signal.time, pysra_sig.values, c=cbox(1))

    sps[1].plot(acc_signal.fa_frequencies, abs(acc_signal.fa_spectrum), c='k')
    sps[1].plot(surf_sig.fa_frequencies, abs(surf_sig.fa_spectrum), c=cbox(0))
    sps[1].plot(pysra_sig.fa_frequencies, abs(pysra_sig.fa_spectrum), c=cbox(1))
    sps[1].set_xlim([0, 20])
    h = surf_sig.smooth_fa_spectrum / acc_signal.smooth_fa_spectrum
    sps[2].plot(surf_sig.smooth_fa_frequencies, h, c=cbox(0))
    pysra_h = pysra_sig.smooth_fa_spectrum / acc_signal.smooth_fa_spectrum
    sps[2].plot(pysra_sig.smooth_fa_frequencies, pysra_h, c=cbox(1))
    sps[2].axhline(1, c='k', ls='--')

    plt.show()
    print(outputs)

if __name__ == '__main__':
    run()