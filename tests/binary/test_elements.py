from tests.binary import functions as fn
import openseespy.opensees as opy


def test_force_beam_column():

    opy.wipe()
    opy.model('basic', '-ndm', 2, '-ndf', 3)  # 2 dimensions, 3 dof per node
    integ_tag = 1
    ele_tag = 1
    fn.define_beam_integration(integ_tag=integ_tag)
    fn.define_force_beam_column(ele_tag, integ_tag)
