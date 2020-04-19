import eqsig
from collections import OrderedDict
import numpy as np
import sfsimodels as sm

import o3seespy as o3



def site_response(sp, asig, linear=0):
    """
    Run seismic analysis of a soil profile - example based on:
    http://opensees.berkeley.edu/wiki/index.php/Site_Response_Analysis_of_a_Layered_Soil_Column_(Total_Stress_Analysis)

    :param sp: sfsimodels.SoilProfile object
        A soil profile
    :param asig: eqsig.AccSignal object
        An acceleration signal
    :return:
    """

    osi = o3.OpenSeesInstance(ndm=2, ndf=2, state=3)
    assert isinstance(sp, sm.SoilProfile)
    sp.gen_split(props=['shear_vel', 'unit_mass', 'cohesion', 'phi', 'bulk_mod', 'poissons_ratio', 'strain_peak'])
    thicknesses = sp.split["thickness"]
    n_node_rows = len(thicknesses) + 1
    node_depths = np.cumsum(sp.split["thickness"])
    node_depths = np.insert(node_depths, 0, 0)
    ele_depths = (node_depths[1:] + node_depths[:-1]) / 2
    shear_vels = sp.split["shear_vel"]
    unit_masses = sp.split["unit_mass"] / 1e3
    g_mods = unit_masses * shear_vels ** 2
    poissons_ratio = sp.split['poissons_ratio']
    youngs_mods = 2 * g_mods * (1 - poissons_ratio)
    bulk_mods = youngs_mods / (3 * (1 - 2 * poissons_ratio))
    bulk_mods = sp.split['bulk_mod'] / 1e3

    ref_pressure = 80.0
    cohesions = sp.split['cohesion'] / 1e3
    phis = sp.split['phi']
    strain_peaks = sp.split['strain_peak']
    grav = 9.81
    damping = 0.03
    omega_1 = 2 * np.pi * 0.5
    omega_2 = 2 * np.pi * 10
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
        nd[f"R{i}L"] = o3.node.Node(osi, 0, -node_depths[i])
        nd[f"R{i}R"] = o3.node.Node(osi, ele_width, -node_depths[i])
        # set x and y dofs equal for left and right nodes
        o3.EqualDOF(osi, nd[f"R{i}L"], nd[f"R{i}R"], [o3.cc.X, o3.cc.Y])

    # Fix base nodes
    o3.Fix2DOF(osi, nd[f"R{n_node_rows - 1}L"], o3.cc.FREE, o3.cc.FIXED)
    o3.Fix2DOF(osi, nd[f"R{n_node_rows - 1}R"], o3.cc.FREE, o3.cc.FIXED)

    # Define dashpot nodes
    dashpot_node_l = o3.node.Node(osi, 0, -node_depths[-1])
    dashpot_node_2 = o3.node.Node(osi, 0, -node_depths[-1])
    o3.Fix2DOF(osi, dashpot_node_l,  o3.cc.FIXED, o3.cc.FIXED)
    o3.Fix2DOF(osi, dashpot_node_2, o3.cc.FREE, o3.cc.FIXED)

    # define equal DOF for dashpot and soil base nodes
    o3.EqualDOF(osi, nd[f"R{n_node_rows - 1}L"], nd[f"R{n_node_rows - 1}R"], [o3.cc.X])
    o3.EqualDOF(osi, nd[f"R{n_node_rows - 1}L"], dashpot_node_2, [o3.cc.X])

    # define materials
    ele_thick = 1.0  # m
    soil_mats = []
    strains = np.logspace(-6, -0.5, 16)
    ref_strain = 0.005
    rats = 1. / (1 + (strains / ref_strain) ** 0.91)
    eles = []
    for i in range(len(thicknesses)):
        if not linear:
            mat = o3.nd_material.PressureIndependMultiYield(osi, 2, unit_masses[i], g_mods[i],
                                                         bulk_mods[i], cohesions[i], strain_peaks[i],
                                                         phis[i], d=0.0, n_surf=16,
                                                         strains=strains, ratios=rats)
        else:
            mat = o3.nd_material.ElasticIsotropic(osi, youngs_mods[i], poissons_ratio[i], rho=unit_masses[i])
        soil_mats.append(mat)

        # def element
        nodes = [nd[f"R{i + 1}L"], nd[f"R{i + 1}R"], nd[f"R{i}R"], nd[f"R{i}L"]]
        ele = o3.element.Quad(osi, nodes, ele_thick, o3.cc.PLANE_STRAIN, mat, b2=grav * unit_masses[i])
        eles.append(ele)

    # define material and element for viscous dampers
    c_base = ele_width * unit_masses[-1] * shear_vels[-1]
    dashpot_mat = o3.uniaxial_material.Viscous(osi, c_base, alpha=1.)
    o3.element.ZeroLength(osi, [dashpot_node_l, dashpot_node_2], mats=[dashpot_mat], dirs=[o3.cc.DOF2D_X])

    # Static analysis
    o3.constraints.Transformation(osi)
    o3.test_check.NormDispIncr(osi, tol=1.0e-4, max_iter=30, p_flag=0)
    o3.algorithm.Newton(osi)
    o3.numberer.RCM(osi)
    o3.system.ProfileSPD(osi)
    o3.integrator.Newmark(osi, newmark_gamma, newmark_beta)
    o3.analysis.Transient(osi)
    o3.analyze(osi, 40, 1.)

    # for i in range(len(soil_mats)):
    #     o3.update_material_stage(osi, soil_mats[i], 1)
    o3.analyze(osi, 50, 0.5)

    # reset time and analysis
    o3.set_time(osi, 0.0)
    o3.wipe_analysis(osi)

    # o3.recorder.NodeToFile(osi, 'sample_out.txt', node=nd["R0L"], dofs=[o3.cc.X], res_type='accel')
    na = o3.recorder.NodeToArrayCache(osi, node=nd["R0L"], dofs=[o3.cc.X], res_type='accel')
    es = o3.recorder.ElementsToArrayCache(osi, eles=eles, arg_vals=['stress'])

    # Define the dynamic analysis
    ts_obj = o3.time_series.Path(osi, dt=asig.dt, values=asig.velocity * -1, factor=c_base)
    o3.pattern.Plain(osi, ts_obj)
    o3.Load(osi, nd["R{0}L".format(n_node_rows - 1)], [1., 0.])

    # Run the dynamic analysis
    o3.algorithm.Newton(osi)
    o3.system.SparseGeneral(osi)
    o3.numberer.RCM(osi)
    o3.constraints.Transformation(osi)
    o3.integrator.Newmark(osi, newmark_gamma, newmark_beta)
    o3.rayleigh.Rayleigh(osi, a0, a1, 0, 0)
    o3.analysis.Transient(osi)

    o3.test_check.EnergyIncr(osi, tol=1.0e-10, max_iter=10)
    analysis_time = asig.time[-1]
    analysis_dt = 0.01
    # o3.extensions.to_py_file(osi)

    while o3.get_time(osi) < analysis_time:
        o3.analyze(osi, 1, analysis_dt)
    o3.wipe(osi)
    outputs = {
        "time": np.arange(0, analysis_time, analysis_dt),
        "rel_disp": [],
        "rel_accel": na.collect(),
        'ele_stresses': es.collect()
    }

    return outputs


def run_pysra(soil_profile, asig, odepths):
    import pysra
    import liquepy as lq
    pysra_profile = lq.sra.sm_profile_to_pysra(soil_profile, d_inc=[0.5] * soil_profile.n_layers)
    # Should be input in g
    pysra_m = pysra.motion.TimeSeriesMotion(asig.label, None, time_step=asig.dt, accels=-asig.values / 9.8)

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


def run(show=0, export=0):
    sl = sm.Soil()

    vs = 160.
    unit_mass = 1700.0
    sl.cohesion = 58.0e3
    sl.phi = 0.0
    sl.g_mod = vs ** 2 * unit_mass
    sl.poissons_ratio = 0.0
    sl.phi = 0.0
    sl.unit_dry_weight = unit_mass * 9.8
    sl.strain_peak = 0.1  # set additional parameter required for PIMY model
    sl.xi = 0.03  # for linear analysis
    assert np.isclose(vs, sl.get_shear_vel(saturated=False))
    soil_profile = sm.SoilProfile()
    soil_profile.add_layer(0, sl)

    sl = sm.Soil()
    vs = 400.
    unit_mass = 1700.0
    sl.g_mod = vs ** 2 * unit_mass
    sl.poissons_ratio = 0.0
    sl.cohesion = 395.0e3
    sl.phi = 0.0
    sl.unit_dry_weight = unit_mass * 9.8
    sl.strain_peak = 0.1  # set additional parameter required for PIMY model
    sl.xi = 0.03  # for linear analysis
    soil_profile.add_layer(9.5, sl)
    soil_profile.height = 20.0
    ecp_out = sm.Output()
    ecp_out.add_to_dict(soil_profile)
    if export:
        import json
        ofile = open('ecp.json', 'w')
        ofile.write(json.dumps(ecp_out.to_dict(), indent=4))
        ofile.close()
    from tests.conftest import TEST_DATA_DIR

    record_path = TEST_DATA_DIR
    record_filename = 'test_motion_dt0p01.txt'
    dt = 0.01
    rec = np.loadtxt(record_path + record_filename) / 2
    acc_signal = eqsig.AccSignal(rec, dt)

    outputs = site_response(soil_profile, acc_signal, linear=0)
    tot_acc = np.sum(abs(outputs['rel_accel']))
    assert np.isclose(tot_acc, 515.76262984), tot_acc  # v3.1.0.11
    resp_dt = outputs['time'][2] - outputs['time'][1]
    surf_sig = eqsig.AccSignal(outputs['rel_accel'], resp_dt)

    if show:
        lw = 0.7
        import matplotlib.pyplot as plt
        from bwplot import cbox
        bf, sps = plt.subplots(nrows=3)

        # linear analysis with pysra
        od = run_pysra(soil_profile, acc_signal, odepths=np.array([0.0, 2.0]))
        pysra_sig = eqsig.AccSignal(od['ACCX_d0'], acc_signal.dt)

        sps[0].plot(acc_signal.time, acc_signal.values, c='k', lw=lw)
        sps[0].plot(surf_sig.time, surf_sig.values, c=cbox(0), lw=lw)
        sps[0].plot(acc_signal.time, pysra_sig.values, c=cbox(1), lw=lw)

        sps[1].plot(acc_signal.fa_frequencies, abs(acc_signal.fa_spectrum), c='k', label='Input', lw=lw)
        sps[1].plot(surf_sig.fa_frequencies, abs(surf_sig.fa_spectrum), c=cbox(0), label='O3', lw=lw)
        sps[1].plot(pysra_sig.fa_frequencies, abs(pysra_sig.fa_spectrum), c=cbox(1), label='pysra', lw=lw)
        sps[1].set_xlim([0, 20])
        h = surf_sig.smooth_fa_spectrum / acc_signal.smooth_fa_spectrum
        sps[2].plot(surf_sig.smooth_fa_frequencies, h, c=cbox(0))
        pysra_h = pysra_sig.smooth_fa_spectrum / acc_signal.smooth_fa_spectrum
        sps[2].plot(pysra_sig.smooth_fa_frequencies, pysra_h, c=cbox(1))
        # o3_nl_h = nl_surf_sig.smooth_fa_spectrum / acc_signal.smooth_fa_spectrum
        # sps[2].plot(nl_surf_sig.smooth_fa_frequencies, o3_nl_h, c=cbox(2))
        sps[2].axhline(1, c='k', ls='--')
        sps[1].legend()
        plt.show()
        print(outputs)


if __name__ == '__main__':
    run(show=1)
