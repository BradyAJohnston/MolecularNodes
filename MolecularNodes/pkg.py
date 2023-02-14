import subprocess
import sys
import os
import site
from importlib.metadata import version
import bpy
import re
import requests
import platform
import pathlib

ADDON_DIR = pathlib.Path(__file__).resolve().parent
print(ADDON_DIR)

PYPI_MIRROR = {
    # the original.
    '':'', 
    # two mirrors in China Mainland to help those poor victims under GFW.
    'BFSU':'https://mirrors.bfsu.edu.cn/pypi/web/simple',
    'TUNA':'https://pypi.tuna.tsinghua.edu.cn/simple',
    # append more if necessary.
}

def verify_user_sitepackages(package_location):
    if os.path.exists(package_location) and package_location not in sys.path:
        sys.path.append(package_location)


def verify(): 
    verify_user_sitepackages(site.getusersitepackages())


def run_pip(cmd, mirror='', timeout=600):
    # path to python.exe
    python_exe = os.path.realpath(sys.executable)
    cmd_list=[python_exe, "-m"] + cmd.split(' ')
    if mirror and mirror.startswith('https'):
        cmd_list+=['-i', mirror]
    try:
        print("Running pip:")
        print(' '.join(cmd_list))
        pip_result = subprocess.run(cmd_list, timeout=timeout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except subprocess.CalledProcessError as e:
        error_message = e.stderr.decode()
        if ("fatal error: 'Python.h' file not found" in error_message) and (platform.system()== "Darwin") and ('arm' in platform.machine()):
            print("BUG: Could not find the Python.h header file in the Blender-build-in-Python.\n" \
                    "This is currently a bug in the Blender of Apple Silicon build.\n" \
                    "Please follow the link to solve it manually: \n" \
                    "https://github.com/BradyAJohnston/MolecularNodes/issues/108#issuecomment-1429384983 ")
        else:
            print("Full error message:\n")
            print(error_message)

def install(pypi_mirror=''):
    # Get PIP upgraded
    run_pip('ensurepip')
    run_pip('pip install --upgrade pip', mirror=PYPI_MIRROR[pypi_mirror])

    #install required packages

    try:
        run_pip(f'pip install -r {ADDON_DIR}/requirements.txt', mirror=PYPI_MIRROR[pypi_mirror])
    except:
        run_pip(f'pip install -r {ADDON_DIR}/requirements.txt', mirror=PYPI_MIRROR['BFSU'])
        

def available():
    verify()
    all_packages_available = True
    for module in ['biotite', 'MDAnalysis']:
        try:
            version(module)
        except Exception as e:
            all_packages_available = False
    return all_packages_available



    
