from logging import error
import sys
import os
import site
import importlib.util

# refreshes the MDAnalysis path and checks if it is installed & importable
# path to python.exe
python_exe = os.path.realpath(sys.executable)


def verify_user_sitepackages(package_location):
    if os.path.exists(package_location) and package_location not in sys.path:
        sys.path.append(package_location)


verify_user_sitepackages(site.getusersitepackages())
verify_user_sitepackages(mda_dir_location)

mda_available = importlib.util.find_spec('MDAnalysis')
if not mda_available:
    error('MDAnalysis is not installed, or not detectable.')