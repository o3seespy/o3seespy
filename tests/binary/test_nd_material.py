import openseespy.opensees as opy
import pytest


def test_elastic_isotropic():
    opy.wipe()
    opy.model('basic', '-ndm', 2, '-ndf', 3)
    v_is_int = 1
    v_is_float = 1.
    # with pytest.raises(opy.error):  # TODO: can't find error class
    #     opy.nDMaterial('ElasticIsotropic', 1, 1., v_is_int, 0.0)
    opy.nDMaterial('ElasticIsotropic', 1, 1., v_is_float, 0.0)

if __name__ == '__main__':
    test_elastic_isotropic()