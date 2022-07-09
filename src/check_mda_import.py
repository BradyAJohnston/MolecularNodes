import sys
import os
import site

# path to python.exe
python_exe = os.path.realpath(sys.executable)


def verify_user_sitepackages(package_location):

    if os.path.exists(package_location) and package_location not in sys.path:
        sys.path.append(package_location)


verify_user_sitepackages(site.getusersitepackages())
verify_user_sitepackages(mda_dir_location)

try:
    import MDAnalysis as mda
    mda_available = True
except:
    mda_available = False
    