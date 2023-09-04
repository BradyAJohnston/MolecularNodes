"""
Handling installation of external python packages inside of Blender's bundled python.
"""

import subprocess
import sys
import os
import logging
from pkg_resources import get_distribution, DistributionNotFound
import bpy
import pathlib
import platform



ADDON_DIR = pathlib.Path(__file__).resolve().parent
"""Folder for the addon on the local machine."""

PYPI_MIRROR = {
    # the original.
    'Default':'', 
    # two mirrors in China Mainland to help those poor victims under GFW.
    'BFSU (Beijing)':'https://mirrors.bfsu.edu.cn/pypi/web/simple',
    'TUNA (Beijing)':'https://pypi.tuna.tsinghua.edu.cn/simple',
    # append more if necessary.
}
"""
Possible PyPi mirrors to install from.
"""

def start_logging(logfile_name: str = 'side-packages-install') -> logging.Logger:
    """
    Configure and start logging to a file.

    Parameters
    ----------
    logfile_name : str, optional
        The name of the log file. Defaults to 'side-packages-install'.

    Returns
    -------
    logging.Logger
        A Logger object that can be used to write log messages.

    This function sets up a logging configuration with a specified log file name and logging level.
    The log file will be created in the `ADDON_DIR/logs` directory. If the directory
    does not exist, it will be created. The function returns a Logger object that can be used to
    write log messages.

    """
    # Create the logs directory if it doesn't exist
    logs_dir = os.path.join(os.path.abspath(ADDON_DIR), 'logs')
    os.makedirs(logs_dir, exist_ok=True)

    # Set up logging configuration
    logfile_path = os.path.join(logs_dir, f"{logfile_name}.log")
    logging.basicConfig(filename=logfile_path, level=logging.INFO)

    # Return logger object
    return logging.getLogger()

"""Determine if the current system is running on Apple Silicon.
True if the system is running on Apple Silicon, False otherwise.
"""
_is_apple_silicon = (sys.platform == "darwin") and ('arm' in platform.machine())

def get_pypi_mirror_alias(self, context, edit_text):
    """
    Get the available PyPI mirror aliases.

    Parameters
    ----------
    self : object
        The object instance.
    context : ...
        The context parameter (description missing).
    edit_text : ...
        The edit_text parameter (description missing).

    Returns
    -------
    keys_view
        A view object of the available PyPI mirror aliases.

    """
    return PYPI_MIRROR.keys()

def process_pypi_mirror_to_url(pypi_mirror_provider: str) -> str: 
    """
    Process a PyPI mirror provider and return the corresponding URL.

    Parameters
    ----------
    pypi_mirror_provider : str
        The PyPI mirror provider to process.

    Returns
    -------
    str
        The URL of the PyPI mirror.

    Raises
    ------
    ValueError
        If the provided PyPI mirror provider is invalid.

    """
    if pypi_mirror_provider.startswith('https:'):
        return pypi_mirror_provider
    elif pypi_mirror_provider in PYPI_MIRROR.keys(): 
        return PYPI_MIRROR[pypi_mirror_provider]
    else:
        raise ValueError(f"Invalid PyPI mirror provider: {pypi_mirror_provider}")


def get_pkgs(requirements: str = None) -> dict:
    """
    Read a requirements file and extract package information into a dictionary.

    Parameters
    ----------
    requirements : str, optional
        The path to the requirements file. If not provided, the function looks for a `requirements.txt`
        file in the same directory as the script.

    Returns
    -------
    dict
        A dictionary containing package information. Each element of the dictionary is a dictionary containing the package name, version, and description.

    Example
    -------
    Given the following requirements file:
    ```python
    Flask==1.1.2 # A micro web framework for Python
    pandas==1.2.3 # A fast, powerful, flexible, and easy-to-use data analysis and manipulation tool
    numpy==1.20.1 # Fundamental package for scientific computing
    ```
    The function would return the following dictionary:
    ```python
    [
        {
            "package": "Flask",
            "version": "1.1.2",
            "desc": "A micro web framework for Python"
        },
        {
            "package": "pandas",
            "version": "1.2.3",
            "desc": "A fast, powerful, flexible, and easy-to-use data analysis and manipulation tool"
        },
        {
            "package": "numpy",
            "version": "1.20.1",
            "desc": "Fundamental package for scientific computing"
        }
    ]
    ```
    """
    import pathlib

    if not requirements:
        folder_path = pathlib.Path(__file__).resolve().parent
        requirements = f"{folder_path}/requirements.txt"

    with open(requirements) as f:
        lines = f.read().splitlines()
        pkgs = {}
        for line in lines:
            try:
                pkg, desc = line.split('#')
                pkg_meta = pkg.split('==')
                name = pkg_meta[0].strip()
                pkgs[name] = {
                    "name": name,
                    "version": pkg_meta[1].strip(),
                    "desc": desc.strip()
                }
            except ValueError:
                # Skip line if it doesn't have the expected format
                pass
    return pkgs

def is_current(package: str) -> bool:
    """
    Check if the specified package is the current version.

    Parameters
    ----------
    package : str
        The name of the package to check.

    Returns
    -------
    bool
        True if the package is the current version, False otherwise.

    """
    pkg = get_pkgs().get(package)
    return is_available(pkg.get('name'), pkg.get('version'))

def is_available(package: str, version: str = None) -> bool:
    """
    Check if a given package is available with the specified version.

    Parameters
    ----------
    package : str
        The name of the package to check.
    version : str, optional
        The version of the package to check.

    Returns
    -------
    bool
        True if the package with the specified version is available, False otherwise.

    Examples
    --------
    >>> is_available('numpy', '1.20.1')
    True
    """

    try: 
        available_version = get_distribution(package).version
        return available_version == version
    except DistributionNotFound:
        return False


def run_python(cmd_list: list=None, mirror_url: str='', timeout: int=600):
    """
    Runs pip command using the specified command list and returns the command output.

    Parameters
    ----------
    cmd_list : list, optional
        List of pip commands to be executed. Defaults to None.
    mirror_url : str, optional
        URL of a package repository mirror to be used for the command. Defaults to ''.
    timeout : int, optional
        Time in seconds to wait for the command to complete. Defaults to 600.

    Returns
    -------
    tuple
        A tuple containing the command list, command return code, command standard output,
        and command standard error.

    Example
    -------
    Install numpy using pip and print the command output
    ```python
    cmd_list = ["-m", "pip", "install", "numpy"]
    mirror_url = 'https://pypi.org/simple/'
    cmd_output = run_python(cmd_list, mirror_url=mirror_url, timeout=300)
    print(cmd_output)
    ```

    """

    # path to python.exe
    python_exe = os.path.realpath(sys.executable)

    # build the command list
    cmd_list=[python_exe] + cmd_list

    # add mirror to the command list if it's valid
    if mirror_url and mirror_url.startswith('https'):
        cmd_list+=['-i', mirror_url]
    
    log = start_logging()
    log.info(f"Running Pip: '{cmd_list}'")

    # run the command and capture the output
    result = subprocess.run(cmd_list, timeout=timeout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if result.returncode != 0:
        log.error('Command failed: %s', cmd_list)
        log.error('stdout: %s', result.stdout.decode())
        log.error('stderr: %s', result.stderr.decode())
    else:
        log.info('Command succeeded: %s', cmd_list)
        log.info('stdout: %s', result.stdout.decode())
    # return the command list, return code, stdout, and stderr as a tuple
    return result

def install_package(package: str, pypi_mirror_provider: str = 'Default') -> list:
    """
    Install a Python package and its dependencies using pip.

    Parameters
    ----------
    package : str
        The name of the package to install.
    pypi_mirror_provider : str, optional
        The name/url of the PyPI mirror provider to use. Default is 'Default'.

    Returns
    -------
    list
        A list of tuples containing the command list, return code, stdout, and stderr
        for each pip command run.

    Raises
    ------
    ValueError
        If the package name is not provided.

    Example
    -------
    To install the package 'requests' from the PyPI mirror 'MyMirror', use:
    ```
    install_package('requests', 'MyMirror')
    ```

    """
    if not package:
        raise ValueError("Package name must be provided.")

    print(f"Installing {package}...")

    mirror_url=process_pypi_mirror_to_url(pypi_mirror_provider=pypi_mirror_provider)
    print(f"Using PyPI mirror: {pypi_mirror_provider} {mirror_url}")
    
    run_python(["-m", "ensurepip"]),
    run_python(["-m", "pip", "install", "--upgrade", "pip"], mirror_url=mirror_url)
    result = run_python(["-m", "pip", "install", package], mirror_url=mirror_url)
    
    return result

class InstallationError(Exception):
    """
    Exception raised when there is an error installing a package.

    Attributes
    ----------
    package_name : str
        The name of the package that failed to install.
    error_message : str
        The error message returned by pip.

    """

    def __init__(self, package_name, error_message):
        self.package_name = package_name
        self.error_message = error_message
        super().__init__(f"Failed to install {package_name}: {error_message}")

def install_all_packages(pypi_mirror_provider: str='Default') -> list:
    """
    Install all packages listed in the 'requirements.txt' file.

    Parameters
    ----------
    pypi_mirror_provider : str, optional
        The PyPI mirror to use for package installation. Defaults to 'Default',
        which uses the official PyPI repository.

    Returns
    -------
    list
        A list of tuples containing the installation results for each package.

    Raises
    ------
    InstallationError
        If there is an error during package installation.

    Example
    -------
    To install all packages listed in the 'requirements.txt' file, run the following command:
    ```
    install_all_packages(pypi_mirror_provider='https://pypi.org/simple/')
    ```

    """
    mirror_url=process_pypi_mirror_to_url(pypi_mirror_provider=pypi_mirror_provider)

    pkgs = get_pkgs()
    results = []
    for pkg in pkgs.items():

        try:
            result = install_package(package=f"{pkg.get('name')}=={pkg.get('version')}", 
                                     pypi_mirror_provider=mirror_url)
            results.append(result)
        except InstallationError as e:
            raise InstallationError(f"Error installing package {pkg.get('name')}: {str(e)}")
    return results

class MN_OT_Install_Package(bpy.types.Operator):
    bl_idname = 'mn.install_package'
    bl_label = 'Install Given Python Package'
    bl_options = {'REGISTER', 'INTERNAL'}
    package: bpy.props.StringProperty(
        name = 'Python Package', 
        description = 'Python Package to Install', 
        default = 'biotite'
    )
    version: bpy.props.StringProperty(
        name = 'Python Package', 
        description = 'Python Package to Install', 
        default = '0.36.1'
    )
    
    description: bpy.props.StringProperty(
        name = 'Operator description', 
        default='Install specified python package.'
    )
    
    @classmethod
    def description(cls, context, properties):
        return properties.description
    
    def execute(self, context):
        installable = f"{self.package}=={self.version}"
        result = install_package(package=installable,
                                 pypi_mirror_provider=bpy.context.scene.pypi_mirror_provider)
        if result.returncode == 0 and is_current(self.package):
            self.report(
                {'INFO'}, 
                f"Successfully installed {self.package} v{self.version}"
                )
        else:
            log_dir = os.path.join(os.path.abspath(ADDON_DIR), 'logs')
            self.report(
                {'ERROR'}, 
                f"Error installing package. Please check the log files in '{log_dir}'."
                )
        return {'FINISHED'}