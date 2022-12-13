# import openseespy.opensees as ops
from o3seespy import opy as ops

def test_basic():
    ops.wipe
    assert 1==1
