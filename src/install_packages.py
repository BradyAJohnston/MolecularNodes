import subprocess
import sys
import os
import site
 
# path to python.exe
python_exe = os.path.realpath(sys.executable)
 
# upgrade pip
subprocess.call([python_exe, "-m", "ensurepip"])
subprocess.call([python_exe, "-m", "pip", "install", "--upgrade", "pip"], timeout=600)
 
# install required packages
subprocess.call([python_exe, "-m", "pip", "install", "atomium"], timeout=600)

def verify_user_sitepackages():
    usersitepackagespath = site.getusersitepackages()

    if os.path.exists(usersitepackagespath) and usersitepackagespath not in sys.path:
        sys.path.append(usersitepackagespath)

verify_user_sitepackages()

try:
    import atomium
    atomium.fetch("1bna")
    atomium_install_successful = True
except:
    atomium_install_successful = False