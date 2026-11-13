import o3seespy.tools

from . import cc, exceptions, extensions, mp, mptools, results

# from .command import test_check  # deprecated
from .__about__ import __version__
from .command import (
    algorithm,
    analysis,
    beam_integration,
    constraints,
    element,
    friction_model,
    geom_transf,
    ground_motion,
    integrator,
    layer,
    mesh,
    nd_material,
    node,
    numberer,
    parallelcmds,
    patch,
    pattern,
    rayleigh,
    recorder,
    region,
    section,
    senscmds,
    system,
    test,
    time_series,
    uniaxial_material,
)
from .command import test as test_check
from .command.common import *
from .opensees_instance import OpenSeesInstance

try:
    from custom_openseespy import opensees as ops
    from custom_openseespy import opensees as opy
except ModuleNotFoundError as error:
    if error.name != "custom_openseespy":
        raise
    try:
        import openseespy.opensees as ops
        import openseespy.opensees as opy
    except ModuleNotFoundError as error:
        if error.name not in ("openseespy", "openseespy.opensees"):
            raise
        import opensees.openseespy as ops
        import opensees.openseespy as opy
