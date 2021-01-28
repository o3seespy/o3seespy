from ..node import Node
from ..element.zero_length import ZeroLength
from .beam_column import ElasticBeamColumn2D
from ..geom_transf import Linear2D as geom_trans_Linear2D
from .bearing import FlatSliderBearing2D
from ..uniaxial_material import Elastic
from ... import cc
from .. import common


class BeamOnNonlinearWinklerFoundation(object):
    def __init__(self, top_nodes, bot_nodes, sf_eles, fd_eles, sf_mats=None, sf_horz_eles=None, sf_horz_mats=None):
        self.top_nodes = top_nodes
        self.bot_nodes = bot_nodes
        self.sf_eles = sf_eles
        self.sf_horz_eles = sf_horz_eles
        self.fd_eles = fd_eles
        self.sf_mats = sf_mats
        self.sf_horz_mats = sf_horz_mats


def gen_shallow_foundation_bnwf(osi, bottom_node, top_node, sf_mats, pos, fd_area, fd_e_mod, fd_iz,
                                    r_flag: float = None, orient: list = None, sf_horz_mats=None, sf_frn=None, offset=0.0):
    """
    Generates nodes and elements for a shallow foundation beam on nonlinear foundation

    Alpha Status - Subject to changes and testing

    Note: Only supports 2D

    Parameters
    ----------
    osi
    bottom_node
    top_node
    sf_mats
    pos
    fd_area
    fd_e_mod
    fd_iz
    r_flag
    orient

    Returns
    -------

    """
    # TODO: sf_horz_mats - if 1 then apply to left node, if len(sf_mats) apply to each
    # TODO: sf_frn - soil-foundation friction, if float, then add each sf_mat inside a flat bearing.
    if not hasattr(sf_mats, '__len__'):  # TODO: only supports 2D
        sf_mats = [sf_mats] * len(pos)
    else:
        sf_mats = list(sf_mats)
    if hasattr(sf_horz_mats, '__len__'):
        assert len(sf_mats) == len(sf_horz_mats)
        sf_horz1_mat = None
    else:
        sf_horz1_mat = sf_horz_mats
        sf_horz_mats = [None] * len(sf_mats)
    injected = 0
    pos = list(pos)
    if offset not in pos:
        pos.append(offset)
        pos.sort()
        injected = 1
    ind = pos.index(offset)
    if injected:
        sf_mats.insert(ind, None)  # insert a null material since there should not be a soil spring here
        sf_horz_mats.insert(ind, None)
    top_nodes = []
    bot_nodes = []
    frn_nodes = []  # unused if no friction model
    con_nodes = top_nodes  # symlink to bot_nodes
    y_inc = 0
    if sf_frn:
        con_nodes = frn_nodes
        y_inc = 0.001  # need small vertical offset for friction element
    sf_eles = []
    sf_horz_eles = []
    fd_eles = []
    transf = geom_trans_Linear2D(osi, [])
    x_centre = bottom_node.x
    y_centre = bottom_node.y
    for i, x in enumerate(pos):
        top_nodes.append(Node(osi, x_centre + x, y_centre))  # node that connects to foundation element
        bot_nodes.append(Node(osi, x_centre + x, y_centre - y_inc))  # node that connects to the far field soil (fixed)
        common.Fix3DOF(osi, bot_nodes[i], cc.FIXED, cc.FIXED, cc.FIXED)
        if sf_frn:  # intermediate nodes to connect friction model between top_nodes and frn_nodes
            frn_nodes.append(Node(osi, x_centre + x, y_centre - y_inc))
            k_add_shear = 1e12
            p_mat = Elastic(osi, 1e12)
            ele = FlatSliderBearing2D(osi, [top_nodes[i], frn_nodes[i]], sf_frn, k_add_shear, p_mat=p_mat,
                                                 mz_mat=p_mat)
        if sf_mats[i] is not None:
            sf_eles.append(ZeroLength(osi, [con_nodes[i], bot_nodes[i]], [sf_mats[i]], [cc.DOF2D_Y], r_flag, orient))
        if sf_horz_mats[i] is not None:
            sf_horz_eles.append(ZeroLength(osi, [con_nodes[i], bot_nodes[i]], [sf_horz_mats[i]], [cc.DOF2D_X], r_flag, orient))
        if i != 0:  # define foundation elements between the top_nodes
            fd_eles.append(ElasticBeamColumn2D(osi, [top_nodes[i - 1], top_nodes[i]], fd_area, fd_e_mod, fd_iz, transf))
    common.EqualDOF(osi, top_node, top_nodes[ind], dofs=[cc.DOF2D_X, cc.DOF2D_Y, cc.DOF2D_ROTZ])
    if sf_horz1_mat:
        sf_horz_mats = ZeroLength(osi, [top_nodes[0], bot_nodes[0]], [sf_horz1_mat], [cc.DOF2D_X], r_flag, orient)
    return BeamOnNonlinearWinklerFoundation(top_nodes, bot_nodes, sf_eles, fd_eles, sf_mats, sf_horz_eles, sf_horz_mats)
