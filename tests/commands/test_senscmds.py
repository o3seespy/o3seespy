import o3seespy as o3  # for testing only
import pytest


def test_parameter():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.senscmds.parameter(osi, p_args=1)

@pytest.mark.skip
def test_add_to_parameter():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.senscmds.add_to_parameter(osi, p_args=1)


@pytest.mark.skip
def test_update_parameter():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.senscmds.update_parameter(osi, new_value=1.0)


@pytest.mark.skip
def test_set_parameter():
    osi = o3.OpenSeesInstance(ndm=2)
    eles = [1, 1]
    o3.senscmds.set_parameter(osi, new_value=1.0, eles=eles, start=1, end=1, args=1)


def test_get_param_tags():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.senscmds.get_param_tags(osi)


@pytest.mark.skip
def test_get_param_value():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.senscmds.get_param_value(osi)


@pytest.mark.skip
def test_compute_gradients():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.senscmds.compute_gradients(osi)


@pytest.mark.skip
def test_sensitivity_algorithm():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.senscmds.sensitivity_algorithm(osi)


@pytest.mark.skip
def test_get_sens_node_disp():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.senscmds.get_sens_node_disp(osi, dof=1, param=[])


@pytest.mark.skip
def test_get_sens_node_vel():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.senscmds.get_sens_node_vel(osi, dof=1, param=[])


@pytest.mark.skip
def test_get_sens_node_accel():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.senscmds.get_sens_node_accel(osi, dof=1, param=[])


@pytest.mark.skip
def test_get_sens_lambda():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.senscmds.get_sens_lambda(osi, param=[])


@pytest.mark.skip
def test_get_sens_section_force():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.senscmds.get_sens_section_force(osi, sec_num=1, dof=1, param=[])


@pytest.mark.skip
def test_get_sens_node_pressure():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.senscmds.get_sens_node_pressure(osi, param=[])

