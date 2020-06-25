from o3seespy.command.element.base_element import ElementBase
from o3seespy.command.node import Node
from o3seespy.command.element.zero_length import ZeroLength
from o3seespy.command.element.beam_column import ElasticBeamColumn2D
from o3seespy.command.geom_transf import Linear2D as geom_trans_Linear2D
from o3seespy import cc
from o3seespy.command import common

class BeamOnNonlinearWinklerFoundation(object):
    def __init__(self, top_nodes, bot_nodes, sf_eles, fd_eles, sf_mats=None):
        self.top_nodes = top_nodes
        self.bot_nodes = bot_nodes
        self.sf_eles = sf_eles
        self.fd_eles = fd_eles
        self.sf_mats = sf_mats


def gen_shallow_foundation_bnwf(osi, bottom_node, top_node, sf_mats, pos, fd_area, fd_e_mod, fd_iz, r_flag: float =None, orient: list =None):
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
        injected = 0
        pos = list(pos)
        if top_node.x not in pos:
            pos.append(top_node.x)
            pos.sort()
            injected = 1
        ind = pos.index(top_node.x)
        if injected:
            sf_mats.insert(ind, None)  # insert a null material since there should not be a soil spring here

        top_nodes = []
        bot_nodes = []
        sf_eles = []
        fd_eles = []
        transf = geom_trans_Linear2D(osi, [])
        x_centre = bottom_node.x
        y_centre = bottom_node.y
        for i, x in enumerate(pos):
            top_nodes.append(Node(osi, x_centre + x, y_centre))
            bot_nodes.append(Node(osi, x_centre + x, y_centre))
            common.Fix3DOF(osi, bot_nodes[i], cc.FIXED, cc.FIXED, cc.FIXED)
            if sf_mats[i] is not None:
                sf_eles.append(ZeroLength(osi, [top_nodes[i], bot_nodes[i]], [sf_mats[i]], [cc.DOF2D_Y], r_flag, orient))
            if i != 0:
                fd_eles.append(ElasticBeamColumn2D(osi, [top_nodes[i-1], top_nodes[i]], fd_area, fd_e_mod, fd_iz, transf))
        common.EqualDOF(osi, top_node, top_nodes[ind], dofs=[cc.DOF2D_X, cc.DOF2D_Y, cc.DOF2D_ROTZ])
        
        return BeamOnNonlinearWinklerFoundation(top_nodes, bot_nodes, sf_eles, fd_eles, sf_mats)
