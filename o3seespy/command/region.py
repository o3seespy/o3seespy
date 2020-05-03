from o3seespy.base_model import OpenSeesObject
from o3seespy.command.common import get_node_tags, get_ele_tags


class NodeRegion(OpenSeesObject):
    op_base_type = "region"
    op_type = None

    def __init__(self, osi, nodes, rayleigh=None):
        """
        A region defined by a group of nodes

        Parameters
        ----------
        osi : opensees_pack.opensees_instance.OpenSeesInstance object
            An instance of opensees
        nodes : str or list
            list of Nodes
        rayleigh: dict
            Dictionary of Rayleigh parameters
        """
        self.osi = osi
        osi.n_region += 1
        self._tag = osi.n_region
        self._parameters = [self._tag, "-node"]

        if isinstance(nodes, str) and nodes == 'all':
            self.nodes = 'all'
            self._parameters += get_node_tags(osi)
        else:
            self.nodes = [x.tag for x in nodes]
            self._parameters += self.nodes
        self.rayleigh = rayleigh
        if rayleigh is not None:
            pms = ['alpha_m', 'beta_k', 'beta_k_init', 'beta_k_comm']
            self._parameters.append('-rayleigh')
            for pm in pms:
                if pm in rayleigh:
                    self._parameters.append(rayleigh[pm])
                else:
                    self._parameters.append(0.0)
        self.to_process(osi)


class ElementRegion(OpenSeesObject):
    op_base_type = "region"
    op_type = None

    def __init__(self, osi, eles, rayleigh=None):
        """
        A region defined by a group of elements

        Parameters
        ----------
        osi : opensees_pack.opensees_instance.OpenSeesInstance object
            An instance of opensees
        eles : str or list
            list of elements
        rayleigh: dict
            Dictionary of Rayleigh parameters
        """
        self.osi = osi
        osi.n_region += 1
        self._tag = osi.n_region
        self._parameters = [self._tag, "-ele"]

        if isinstance(eles, str) and eles == 'all':
            self.eles = 'all'
            self._parameters += get_ele_tags(osi)
        else:
            self.eles = [x.tag for x in eles]
            self._parameters += self.eles
        self.rayleigh = rayleigh
        if rayleigh is not None:
            pms = ['alpha_m', 'beta_k', 'beta_k_init', 'beta_k_comm']
            self._parameters.append('-rayleigh')
            for pm in pms:
                if pm in rayleigh:
                    self._parameters.append(rayleigh[pm])
                else:
                    self._parameters.append(0.0)
        self.to_process(osi)
