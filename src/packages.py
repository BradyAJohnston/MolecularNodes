import subprocess
import sys
import os
import site

def verify_user_sitepackages(package_location):
    if os.path.exists(package_location) and package_location not in sys.path:
        sys.path.append(package_location)


def verify(): 
    verify_user_sitepackages(site.getusersitepackages())

def install_required_packages():
    # path to python.exe
    python_exe = os.path.realpath(sys.executable)

    # upgrade pip
    subprocess.call([python_exe, "-m", "ensurepip"])
    subprocess.call([python_exe, "-m", "pip", "install","--upgrade", "pip"], timeout=600)

    #install required packages
    subprocess.call([python_exe, "-m", "pip", "install", "biotite"], timeout=600)

success = False
def install_packages():
    try:
        import biotite
        success = True
    except: 
        verify()
    try:
        import biotite
        success = True
    except:
        install_required_packages()
    try: 
        import biotite
        success = True
    except: 
        print("Still Can't import Biotite")