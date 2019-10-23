from openseespy import opensees as op


def elastic_section(section_id):

    e_conc = 30e6
    area = 0.3 * 0.4
    inertia = 0.3 * 0.4 ** 3 / 12

    op.section("Elastic", section_id, *[e_conc, area, inertia])


def bilinear_material(mat_id):
    e_conc = 30e6
    depth = 0.4
    width = 0.3
    inertia = width * depth ** 3 / 12
    ei = e_conc * inertia
    eps_yield = 300.0e6 / 200e9
    phi_y = 2.1 * eps_yield / depth
    mat_props = [ei, 0.05 * ei, phi_y]
    op.uniaxialMaterial("ElasticBilin", mat_id, *mat_props)


def uniaxial_steel01_material(mat_id):
    mat_type = 'Steel01'

    mat_args = [300e6, 200e9, 0.001]
    op.uniaxialMaterial(mat_type, mat_id, *mat_args)


def uniaxial_steel01_section(section_id, mat_id=None):
    if mat_id is None:
        mat_id = section_id
    uniaxial_steel01_material(mat_id)

    op.section("Uniaxial", section_id, mat_id, "Mz")


def define_beam_integration(integ_tag):
    lp_i = 0.1
    lp_j = 0.2
    section1_id = 1
    section2_id = 2
    section3_id = 3
    uniaxial_steel01_section(section1_id)
    uniaxial_steel01_section(section2_id)
    uniaxial_steel01_section(section3_id)
    op.beamIntegration('HingeMidpoint', integ_tag, section1_id, lp_i, section2_id, lp_j, section3_id)


def define_force_beam_column(ele_i, integ_tag):

    # Establish nodes
    bot_node = 1
    top_node = 2
    transf_tag = 1
    op.geomTransf('Linear', transf_tag, *[])
    op.node(bot_node, 0., 0.)
    op.node(top_node, 0., 5.)
    op.element('forceBeamColumn', ele_i, bot_node, top_node, transf_tag, integ_tag)
