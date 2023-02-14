import subprocess
import sys
import os
import site
from importlib.metadata import version

packages = [
    "biotite==0.36.1", 
    "MDAnalysis==2.2.0", 
    "starfile"
    ]

def verify_user_sitepackages(package_location):
    if os.path.exists(package_location) and package_location not in sys.path:
        sys.path.append(package_location)


def verify(): 
    verify_user_sitepackages(site.getusersitepackages())

def install():
    # path to python.exe
    python_exe = os.path.realpath(sys.executable)
    # upgrade pip
    subprocess.call([python_exe, "-m", "ensurepip"])
    subprocess.call([python_exe, "-m", "pip", "install","--upgrade", "pip"], timeout=600)
    #install required packages
    for pkg in packages:
        subprocess.call([python_exe, "-m", "pip", "install", pkg], timeout=600)

def available():
    verify()
    all_packages_available = True
    for module in ['biotite', 'MDAnalysis', 'starfile']:
        try:
            version(module)
        except:
            all_packages_available = False
    return all_packages_available