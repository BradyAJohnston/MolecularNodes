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

# for M1 installation issue. 
# the keys should be exactly as those in `PYPI_MIRROR`
CONDA_MIRROR_MAC_ARM = {
    # the original.
    'Default': 'https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh',
    # two mirrors in China Mainland to help those poor victims under GFW.
    'BFSU (Beijing)': 'https://mirrors.bfsu.edu.cn/anaconda/miniconda/Miniconda3-latest-MacOSX-arm64.sh',
    'TUNA (Beijing)': 'https://mirrors.bfsu.edu.cn/anaconda/miniconda/Miniconda3-latest-MacOSX-arm64.sh'
    # append more if necessary.
}

def does_it_ARMed(binary_path):
    assert pathlib.Path.exists(binary_path)
    test_cmd_run=subprocess.run(['file',binary_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout=test_cmd_run.stdout.decode().strip()
    stderr=test_cmd_run.stderr.decode().strip()
    if test_cmd_run.returncode != 0:
        return ('arm64' in stdout)
    else:
        raise RuntimeError(f'Failed on determine the architechture of {binary_path}: \nFull error: {stderr}')
    
        

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
        pip_result = subprocess.run(cmd_list, timeout=timeout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,env=env)
        return [cmd_list,pip_result.stdout.decode()]    
    except subprocess.CalledProcessError as e:
        error_message = e.stderr.decode()
        return [cmd_list,"Full error message:\n" + error_message]

def install(mirror=''):
    import platform
    install_commands=[]
    install_logs=[]

    CPATH_ANOUNCEMENT={}

    if (platform.system()== "Darwin") and ('arm' in platform.machine()):
        # M1 issue
        # find default conda path and py310
        
        TEST_SYSTEM_PYTHON310_CFG = subprocess.run(['which','python3.10-config'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        SYSTEM_PYTHON310=pathlib.Path(TEST_SYSTEM_PYTHON310_CFG.stdout.decode().strip().replace('-config','')).resolve() if TEST_SYSTEM_PYTHON310_CFG.returncode==0 else ''

        # basic pkg managements
        CONDA_PATH_EXEC=pathlib.Path(os.environ['CONDA_EXE']).resolve() if 'CONDA_EXE' in os.environ.keys() else ''
        install_commands.append('CONDA detection..')
        install_logs.append(CONDA_PATH_EXEC)
        
        
        if not (SYSTEM_PYTHON310 and does_it_ARMed(SYSTEM_PYTHON310)):
            if not CONDA_PATH_EXEC:
                print(f'Conda is required for dependency installation, as the \'Python.h\' header file is missing in Python bundled with the Apple Silicon versions of Blender.')

                # install a new miniconda and apply it
                MINICONDA_INSTALLER_URL = CONDA_MIRROR_MAC_ARM[mirror]
                MINICONDA_INSTALLER_PATH = '/tmp/miniconda_installer.sh'
                MINICONDA_INSTALL_PATH=pathlib.Path(sys.executable).resolve().parent.parent.parent.joinpath('miniconda_mn_py310')
                
                # download miniconda
                miniconda_download_cmd=['curl', '-o', MINICONDA_INSTALLER_PATH, '-L', MINICONDA_INSTALLER_URL]
                install_commands.append(miniconda_download_cmd)
                download_miniconda_run=subprocess.run(miniconda_download_cmd, check=True,stdout=subprocess.PIPE)
                install_logs.append(f'Miniconda is dowloaded from {MINICONDA_INSTALLER_URL} and will be stored at {MINICONDA_INSTALLER_PATH}\n{download_miniconda_run.stdout.decode()}\n')
                
                # install miniconda
                miniconda_install_cmd=['sh', MINICONDA_INSTALLER_PATH, '-b','-f', '-p', MINICONDA_INSTALL_PATH]
                install_commands.append(miniconda_install_cmd)
                install_miniconda_run=subprocess.run(miniconda_install_cmd, check=True,stdout=subprocess.PIPE)
                CONDA_PATH_EXEC = pathlib.Path(f'{MINICONDA_INSTALL_PATH}/bin/conda').resolve()
                install_logs.append(f'Miniconda has been installed successfully to {MINICONDA_INSTALL_PATH}.\n{install_miniconda_run.stdout.decode()}\n')

                # apply new conda exe to PATH
                os.environ['CONDA_EXE'] = str(CONDA_PATH_EXEC)
                os.environ['PATH'] = f'{str(CONDA_PATH_EXEC.parent)}:{os.environ["PATH"]}'
                install_commands.append('Set up environment variables..')
                install_logs.append(f'CONDA_EXE={os.environ["CONDA_EXE"]}, PATH={os.environ["PATH"]}')
        

            CONDA_PATH=CONDA_PATH_EXEC.parent.parent
            CONDA_MN_PREFIX='MN_PY310'
            SUPPOSED_PYTHONPATH=CONDA_PATH.joinpath('envs',CONDA_MN_PREFIX,'bin','python3.10')
            install_commands.append(f'Detecting Conda env with python3.10..')

            # if a conda env has already been created before.
            if pathlib.Path.exists(SUPPOSED_PYTHONPATH):
                install_logs.append(SUPPOSED_PYTHONPATH)
                print(f'Found previously installed Python 3.10 via Miniconda: {SUPPOSED_PYTHONPATH}')
                SYSTEM_PYTHON310=SUPPOSED_PYTHONPATH
            else:
                install_logs.append('Not found...')
                # run a new conda environment setup
                conda_cmd_new_env=[ CONDA_PATH_EXEC,'create','-n',CONDA_MN_PREFIX, 'python=3.10', '-y']
                conda_install_run=subprocess.run(conda_cmd_new_env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env={'CONDA_SUBDIR':'osx-arm64'})
                install_commands.append(conda_cmd_new_env)
                install_logs.append(conda_install_run.stdout.decode() + conda_install_run.stderr.decode())
                
                # check the conda installation
                if conda_install_run.returncode == 0 and pathlib.Path.exists(SUPPOSED_PYTHONPATH):
                    SYSTEM_PYTHON310=SUPPOSED_PYTHONPATH
                    print(f'Python 3.10 installed: {SYSTEM_PYTHON310}')
                else:
                    print(f'Failed to setup new conda environment {CONDA_MN_PREFIX}.')
                    print(conda_install_run.stderr.decode())
            
            # now we determine the python cpath include
            if SYSTEM_PYTHON310:
                system_python310_cmd=[f'{SYSTEM_PYTHON310}-config','--include']
                SYSTEM_PYTHON310_INCLUDED=subprocess.run(system_python310_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,).stdout.decode().strip()
                install_commands.append(system_python310_cmd)
                install_logs.append(SYSTEM_PYTHON310_INCLUDED)
                # reading the results
                CPATH_ANOUNCEMENT={'CPATH':SYSTEM_PYTHON310_INCLUDED.replace('-I',':').replace(' ',''),'PATH':os.environ['PATH']}
                print(f'PIP w/ ENV: {CPATH_ANOUNCEMENT}')
                install_commands.append("CPATH_ANOUNCEMENT")
                install_logs.append(CPATH_ANOUNCEMENT)
            else:
                print(f'Python 3.10 is not installed before dependency installation.')
    

    for cmd_list,stdouterr in [
        # Get PIP upgraded
        run_pip('ensurepip'),
        run_pip('pip install --upgrade pip', mirror=PYPI_MIRROR[mirror]),
        # dependencies
        run_pip(['pip', 'install', '-r', f'{ADDON_DIR}/requirements.txt'], mirror=PYPI_MIRROR[mirror],env=CPATH_ANOUNCEMENT)]:
        
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
                logfile.write(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n{cmd}\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n{log}\n\n')

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
                message=f'Failed to install required packages. \nPlease check log file for details: {logfile_path}'
                )
        
        return {'FINISHED'}