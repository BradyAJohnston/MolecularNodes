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


def run_pip(cmd, mirror='', timeout=600):
    # path to python.exe
    python_exe = os.path.realpath(sys.executable)
    if type(cmd)==list:
        cmd_list=[python_exe, "-m"] + cmd
    elif type(cmd)==str:
        cmd_list=[python_exe, "-m"] + cmd.split(" ")
    else:
        raise TypeError(f"Invalid type of input cmd.")
    if mirror and mirror.startswith('https'):
        cmd_list+=['-i', mirror]
    try:
        print("Running pip:")
        print(cmd_list)
        pip_result = subprocess.run(cmd_list, timeout=timeout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return [cmd_list,pip_result.stdout.decode()]    
    except subprocess.CalledProcessError as e:
        error_message = e.stderr.decode()
        return [cmd_list,"Full error message:\n" + error_message]

def install(pypi_mirror=''):
    
    install_commands=[]
    install_logs=[]

    # Get PIP upgraded


    for cmd_list,stdouterr in [
        run_pip('ensurepip'),
        run_pip('pip install --upgrade pip', mirror=PYPI_MIRROR[pypi_mirror]),
        run_pip(cmd=['pip', 'install', '-r', f'{ADDON_DIR}/requirements.txt'], mirror=PYPI_MIRROR[pypi_mirror])]:
        
        install_commands.append(cmd_list)
        install_logs.append(stdouterr)
    return [install_commands,install_logs]

def available():
    verify()
    all_packages_available = True
    for module in ['biotite', 'MDAnalysis']:
        try:
            version(module)
        except Exception as e:
            all_packages_available = False
    return all_packages_available

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
            install_commands,install_logs=install(
                pypi_mirror=bpy.context.scene.pypi_mirror,
            )

            # log cmd and stdout/stderr
            for cmd,log in zip(install_commands,install_logs):
                logfile.write(f'{cmd}\n{log}\n')

                # solve this Apple Silicon installation failure 
                if ("fatal error: 'Python.h' file not found" in log) and (platform.system()== "Darwin") and ('arm' in platform.machine()):
                    error_msg = f"ERROR: Could not find the 'Python.h' header file in version of Python bundled with Blender.\n" \
                            "This is a problem with the Apple Silicon versions of Blender.\n" \
                            "Please follow the link to the MolecularNodes GitHub page to solve it manually: \n" \
                            "https://github.com/BradyAJohnston/MolecularNodes/issues/108#issuecomment-1429384983 "
                    logfile.write(f"{error_msg}\n")
                    # print this message to Blender script window
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
            
        else:
            # bpy.context.preferences.addons['MolecularNodesPref'].preferences.packages_available = False
            self.report(
                {'ERROR'}, 
                message=f'Failed to install required packages. \n {error_msg}. \nPlease check log file for details: {logfile_path}'
                )
        
        return {'FINISHED'}