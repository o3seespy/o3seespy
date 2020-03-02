import sys
from tests import conftest

sys.path.append(conftest.EXAMPLES_DIR)
import cantilever_w_hinge_under_force_then_disp_loading
import frame_building_dynamic_analysis
import readme_example
import simple_sdof
import site_response_analysis
import pytest


def test_cantilever_w_hinge_under_force_then_disp_loading():
    cantilever_w_hinge_under_force_then_disp_loading.run(use_pload=0)


def test_frame_building_dynamic_analysis():
    frame_building_dynamic_analysis.run_as_e2e()


def test_readme_example():
    readme_example.run(show=0)


def test_simple_sdof():
    simple_sdof.test_sdof(show=0)


def test_site_response_analysis():
    site_response_analysis.run(show=0)
