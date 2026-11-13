import subprocess
import sys
import unittest
from pathlib import Path

IMPORT_CHECK = """
import builtins
import sys
import types

selected, missing_dependency = sys.argv[1:]
packages = ('custom_openseespy', 'openseespy', 'opensees')
backends = {package: types.ModuleType(package + '_backend') for package in packages}
backend = backends[selected]
real_import = builtins.__import__

class BackendError(Exception):
    pass

for package in ('custom_openseespy', 'openseespy'):
    backends[package].OpenSeesError = BackendError

def import_backend(name, globals=None, locals=None, fromlist=(), level=0):
    if level == 0 and name == 'opensees.errors':
        return types.SimpleNamespace(XaraError=BackendError)
    if level == 0 and name in ('custom_openseespy', 'openseespy.opensees', 'opensees.openseespy'):
        package = name.split('.')[0]
        if package == missing_dependency:
            raise ModuleNotFoundError('missing internal dependency', name='internal_dependency')
        if packages.index(package) < packages.index(selected):
            raise ModuleNotFoundError('missing backend', name=package)
        module = types.ModuleType(package)
        setattr(module, 'openseespy' if package == 'opensees' else 'opensees', backends[package])
        return module
    return real_import(name, globals, locals, fromlist, level)

builtins.__import__ = import_backend
try:
    import o3seespy
except ModuleNotFoundError as error:
    assert missing_dependency and error.name == 'internal_dependency', error
else:
    assert not missing_dependency, 'expected a missing internal dependency error'
    import o3seespy.base_model
    import o3seespy.opensees_instance
    for module in (o3seespy, o3seespy.base_model, o3seespy.opensees_instance):
        assert module.opy is backend and module.ops is backend, module.__name__

    def fail(*args):
        raise BackendError('backend failure')

    backend.node = fail
    instance = o3seespy.OpenSeesInstance.__new__(o3seespy.OpenSeesInstance)
    obj = o3seespy.base_model.OpenSeesObject()
    obj.op_base_type = 'node'
    obj._parameters = [1]
    for call in (lambda: instance.to_opensees('node', [1]), obj.to_opensees):
        try:
            call()
        except ValueError as error:
            assert 'backend failure' in str(error), error
        else:
            raise AssertionError('expected ValueError')
"""


class TestBackendImports(unittest.TestCase):
    def test_backend_selection_and_internal_errors(self):
        cases = (
            ("custom_openseespy", ""),
            ("openseespy", ""),
            ("opensees", ""),
            ("custom_openseespy", "custom_openseespy"),
            ("openseespy", "openseespy"),
        )
        for selected, missing_dependency in cases:
            with self.subTest(selected=selected, missing_dependency=missing_dependency):
                result = subprocess.run(
                    [sys.executable, "-c", IMPORT_CHECK, selected, missing_dependency],
                    cwd=Path(__file__).resolve().parents[1],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
