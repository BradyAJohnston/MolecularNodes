ADDON = "MolecularNodes"

import os
import sys
import MolecularNodes as mn

mn.pkg.install_package('pytest-snapshot')

try:
    import pytest
except Exception as e:
    print(e)
    sys.exit(1)

try:
    sys.path.append(os.environ["LOCAL_PYTHONPATH"])
    from addon_helper import zip_addon, change_addon_dir, install_addon, cleanup
except Exception as e:
    print(e)
    sys.exit(1)


class SetupPlugin:
    def __init__(self, addon):
        self.addon = addon
        self.addon_dir = "MolecularNodes"

    def pytest_configure(self, config):
        (self.bpy_module, self.zfile) = zip_addon(self.addon, self.addon_dir)
        change_addon_dir(self.bpy_module, self.addon_dir)
        install_addon(self.bpy_module, self.zfile)
        config.cache.set("bpy_module", self.bpy_module)

    def pytest_unconfigure(self):
#         cmd = "coverage xml"
#         os.system(cmd)
        cleanup(self.addon, self.bpy_module, self.addon_dir)
        print("*** test run reporting finished")


try:
    exit_val = pytest.main(["tests", "-v", "-x", "--cov", "--cov-report", "term-missing", "--cov-report", "xml",], plugins=[SetupPlugin(ADDON)])
except Exception as e:
    print(e)
    exit_val = 1
sys.exit(exit_val)
