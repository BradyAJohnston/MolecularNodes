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



def get_mirror_alias(self, context, edit_text):
    return PYPI_MIRROR.keys()
    

def verify_user_sitepackages(package_location):
    if os.path.exists(package_location) and package_location not in sys.path:
        sys.path.append(package_location)


def verify(): 
    verify_user_sitepackages(site.getusersitepackages())


def run_pip(cmd, mirror='', timeout=600,env={}):
    # path to python.exe
    python_exe = os.path.realpath(sys.executable)
    assert type(cmd)==list # TODO remove this assertion via typing check in future refactoring
    cmd_list=[python_exe, "-m"] + cmd
    if mirror and mirror.startswith('https'):
        cmd_list+=['-i', mirror]
    try:
        print("Running pip:")
        print(cmd_list)
        pip_result = subprocess.run(cmd_list, timeout=timeout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,env=env)
        return [cmd_list,pip_result.stdout.decode()]    
    except subprocess.CalledProcessError as e:
        error_message = e.stderr.decode()
        return [cmd_list,"Full error message:\n" + error_message]

def install(mirror=''):
    import platform
    install_commands=[]
    install_logs=[]

    # check if a given binary is built in ARM64 arch for MacOS
    # This will be used only when M1/M2 issue occurs.
    def does_it_ARMed(binary_path):
        assert pathlib.Path.exists(binary_path)
        test_cmd_run=subprocess.run(['file',binary_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout=test_cmd_run.stdout.decode().strip()
        stderr=test_cmd_run.stderr.decode().strip()
        if test_cmd_run.returncode == 0:
            return ('arm64' in stdout)
        else:
            print(f'Failed on determine the architechture of {binary_path}: \nFull error: {stderr}')
            return False


    CPATH_ANOUNCEMENT={}

    # M1 issue caused by missing 'Python.h'(python-dev) in bundled python package of Blender itself.
    # to solve this problem, one posible solution is to borrow the header file from external python3.10 CPATH
    # please see the phenotype solution at https://github.com/BradyAJohnston/MolecularNodes/issues/108#issuecomment-1429384983 
    # and a shell helper script in https://github.com/BradyAJohnston/MolecularNodes/issues/108#issuecomment-1442794000.
    # in this circunstance, system-wide Python is firstly checked by looking for `python3.10-config`. 
    # If it doesn't exist, an error will be raised.
    # Since GCC is required for building MDAnaysis(Numpy?) from scratch, Xcode commandline tools should be installed before this installation.
    if (platform.system()== "Darwin") and ('arm' in platform.machine()):
        # good luck guys.

        # detect py310
        TEST_SYSTEM_PYTHON310_CFG = subprocess.run(['which','python3.10-config'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        SYSTEM_PYTHON310=pathlib.Path(TEST_SYSTEM_PYTHON310_CFG.stdout.decode().strip().replace('-config','')).resolve() if TEST_SYSTEM_PYTHON310_CFG.returncode==0 else ''
        install_commands.append('System-wide detection..')
        install_logs.append(f"{(SYSTEM_PYTHON310 != '')}: {SYSTEM_PYTHON310}")

        # now we determine the python cpath include
        if SYSTEM_PYTHON310 and does_it_ARMed(SYSTEM_PYTHON310):
            system_python310_cmd=[f'{SYSTEM_PYTHON310}-config','--include']
            install_commands.append(f'Detect External Python include path: {system_python310_cmd}')
            SYSTEM_PYTHON310_INCLUDED=subprocess.run(system_python310_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,).stdout.decode().strip()
            install_logs.append(SYSTEM_PYTHON310_INCLUDED)

            # PATH is required by biopython for (--use-pep517) issue (?)
            # gcc is required for building the dependency MDAnaysis from scratch.
            CPATH_ANOUNCEMENT={'CPATH':SYSTEM_PYTHON310_INCLUDED.replace('-I',':').replace(' ',''),'PATH':os.environ['PATH']}
            install_commands.append("Final CPATH_ANOUNCEMENT for M1/M2 issue")
            install_logs.append(CPATH_ANOUNCEMENT)
        else:
            # TODO: add guideline url to both logfile and error message.
            install_commands.append("Error: ")
            install_logs.append("System-wide Python3.10 should be available. Please see the GitHub page of MolecularNodes for more information.")
            # If system-wide Python3.10 in arm64 is not available, raise an error.
            raise RuntimeError("System-wide Python3.10 should be available. Please see the GitHub page of MolecularNodes for more information.")
    
    for cmd_list,stdouterr in [
        # Get PIP installed and upgraded
        run_pip(['ensurepip']),
        run_pip(['pip', 'install', '--upgrade','pip'], mirror=PYPI_MIRROR[mirror]),
        # dependencies
        run_pip(['pip', 'install', '-r', f'{ADDON_DIR}/requirements.txt'], mirror=PYPI_MIRROR[mirror], env=CPATH_ANOUNCEMENT)]:
        
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
            import datetime
            
            # generate logfile
            logfile_path = os.path.abspath(str(ADDON_DIR) + "/logs/" + 'side-packages-install.log')
            logfile = open(logfile_path, 'a')
            
            logfile.write("-----------------------------------" + '\n')
            logfile.write("Installer Started: " + str(datetime.datetime.now()) + '\n')
            logfile.write("-----------------------------------" + '\n')
            install_commands,install_logs=install(
                mirror=bpy.context.scene.mirror,
            )

            # log cmd and stdout/stderr
            for cmd,log in zip(install_commands,install_logs):
                logfile.write(
                            '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n' \
                            f'{cmd}\n' \
                            '===================================\n' \
                            f'{log}\n'\
                            '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n')

            logfile.write("###################################" + '\n')
            logfile.write("Installer finished: " + str(datetime.datetime.now()) + '\n')
            logfile.write("###################################" + '\n')
            
            # close the logfile
            logfile.close()

        if available():
            # bpy.context.preferences.addons['MolecularNodesPref'].preferences.packages_available = True
            self.report(
                {'INFO'}, 
                message='Successfully Installed Required Packages.'
                )
            
        else:
            # bpy.context.preferences.addons['MolecularNodesPref'].preferences.packages_available = False
            self.report(
                {'ERROR'}, 
                message=f'Failed to install required packages. \nPlease check log file for details: {logfile_path} .'
                )
        
        return {'FINISHED'}