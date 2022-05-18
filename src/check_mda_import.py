import sys
import os
import site

# path to python.exe
python_exe = os.path.realpath(sys.executable)


def verify_user_sitepackages(usersitepackagespath):
    usersitepackagespath = site.getusersitepackages()

    if os.path.exists(usersitepackagespath) and usersitepackagespath not in sys.path:
        sys.path.append(usersitepackagespath)


verify_user_sitepackages(site.getusersitepackages())
verify_user_sitepackages(mdanalysis_dir_location)

try:
    import MDAnalysis as mda
    mda_available = True
except:
    mda_available = False
