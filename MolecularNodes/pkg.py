import subprocess
import sys
import os
import site
from importlib.metadata import version
import bpy
import pathlib

ADDON_DIR = pathlib.Path(__file__).resolve().parent
print(ADDON_DIR)

PYPI_MIRROR = {
    # the original.
    'Default':'', 
    # two mirrors in China Mainland to help those poor victims under GFW.
    'BFSU (Beijing)':'https://mirrors.bfsu.edu.cn/pypi/web/simple',
    'TUNA (Beijing)':'https://pypi.tuna.tsinghua.edu.cn/simple',
    # append more if necessary.

}

def get_pypi_mirror_alias(self, context, edit_text):
    return PYPI_MIRROR.keys()
    

def verify_user_sitepackages(package_location):
    if os.path.exists(package_location) and package_location not in sys.path:
        sys.path.append(package_location)


def verify(): 
    verify_user_sitepackages(site.getusersitepackages())


def run_pip(cmd_list=[], mirror='', timeout=600):
    # path to python.exe
    python_exe = os.path.realpath(sys.executable)

    cmd_list=[python_exe, "-m"] + cmd_list

    if mirror and mirror.startswith('https'):
        cmd_list+=['-i', mirror]

    print("Running pip:")
    print(cmd_list)
    pip_result = subprocess.run(cmd_list, timeout=timeout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return (cmd_list,pip_result.returncode,pip_result.stdout.decode(),pip_result.stderr.decode())  

def install(pypi_mirror=''):

    results=[]


    for result in [
        # Get PIP upgraded
        run_pip(['ensurepip']),
        run_pip(['pip', 'install', '--upgrade','pip'], mirror=PYPI_MIRROR[pypi_mirror]),
        # Install from requirements
        run_pip(['pip', 'install', '-r', f'{ADDON_DIR}/requirements.txt'], mirror=PYPI_MIRROR[pypi_mirror])]:
        results.append(result)
    return results


def available():
    import pkg_resources
    verify()
    # compare the pinned versions of requirements.txt against the installed.
    with open(f'{ADDON_DIR}/requirements.txt') as f:
        requirements = f.read().splitlines()

    for requirement in requirements:
        if requirement.startswith('#') or requirement.strip() == '': continue
        req = pkg_resources.Requirement.parse(requirement)
        try:
            version(req.name)
            installed = pkg_resources.get_distribution(req.name)
            # print to VSC console for debug
            print(f'Required: {req}\tInstalled: {installed}\tMatch: {installed.version in req}')
            if not installed.version in req:
                return False
        except Exception as e:
            print(e)
            return False

    return True


class MOL_OT_install_dependencies(bpy.types.Operator):
    bl_idname = "mol.install_dependencies"
    bl_label = "Install Dependencies"
    bl_description = "Install the required python packages to enable import."
    bl_options = {'REGISTER', 'INTERNAL'}
    
    def execute(self, context):
        if not available():
            import datetime,platform
            
            # generate logfile
            logfile_path = os.path.abspath(str(ADDON_DIR) + "/logs/" + 'side-packages-install.log')
            logfile = open(logfile_path, 'a')
            
            logfile.write("-----------------------------------" + '\n')
            logfile.write("Installer Started: " + str(datetime.datetime.now()) + '\n')
            logfile.write("-----------------------------------" + '\n')
            install_results=install(
                pypi_mirror=bpy.context.scene.pypi_mirror,
            )
            
            # a lable to check all return codes of installation results.
            no_errors=True

            # log cmd, return code, and stdout/stderr
            for (cmd, returncode, stdout, stderr) in install_results:
                logfile.write(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n{cmd}\n{returncode}\n===================================\n{stdout}\n===================================\n{stderr}\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n')
                error_msg=''

                # for dependency updates. if all return code is 0, then update is defined as successful.
                # the furture call of available() will still report False, unless the Blender is restarted to reload MN
                # no_errors is used here to check if there are any errors in the latest pip installation.
                if returncode != 0: no_errors=False

                # provide current solution for this Apple Silicon installation failure 
                if ("fatal error: 'Python.h' file not found" in stderr) and (platform.system()== "Darwin") and ('arm' in platform.machine()):
                    error_msg = f"ERROR: Could not find the 'Python.h' header file in version of Python bundled with Blender.\n" \
                            "This is a problem with the Apple Silicon versions of Blender.\n" \
                            "Please follow the link to the MolecularNodes GitHub page to solve it manually: \n" \
                            "https://github.com/BradyAJohnston/MolecularNodes/issues/108#issuecomment-1429384983 "
                    logfile.write(f"{error_msg}\n")

                    # print this message to VSC console window
                    print(error_msg)


            logfile.write("###################################" + '\n')
            logfile.write("Installer finished: " + str(datetime.datetime.now()) + '\n')
            logfile.write("###################################" + '\n')
            
            # close the logfile
            logfile.close()
        
        if available():
            # bpy.context.preferences.addons['MolecularNodesPref'].preferences.packages_available = True
            self.report(
                {'INFO'}, 
                message='Successfully Installed Required Packages'
                )
            
        elif no_errors:
            # update via requirements. Ask user to restart Blender for a complete reload.
            self.report(
                {'ERROR'}, 
                message=f'Successfully Installed Required Packages.\nPlease restart Blender and check log file for details: {logfile_path}'
                )
            
        else:
            # bpy.context.preferences.addons['MolecularNodesPref'].preferences.packages_available = False
            self.report(
                {'ERROR'}, 
                message=f'Failed to install required packages. \n{error_msg}. \nPlease check log file for details: {logfile_path}'
                )
        
        return {'FINISHED'}