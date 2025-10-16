import subprocess
import sys


def install_dependency_group(pyproject_path, group_name):
    """
    Install packages from a specific dependency group in pyproject.toml.

    Args:
        pyproject_path: Path to pyproject.toml file
        group_name: Name of the dependency group to install
    """
    python_exe = sys.executable
    subprocess.check_call([python_exe, "-m", "pip", "install", "uv"])
    subprocess.check_call([python_exe, "-m", "pip", "install", "-e", "."])
    subprocess.check_call(
        [
            python_exe,
            "-m",
            "uv",
            "pip",
            "install",
            "-r",
            pyproject_path,
            "--extra",
            group_name,
        ]
    )


# Configuration
PYPROJECT_PATH = "pyproject.toml"  # Change this to your file path
GROUP_NAME = "test"  # Change this to your group name (e.g., "dev", "test", "docs")

if __name__ == "__main__":
    install_dependency_group(PYPROJECT_PATH, GROUP_NAME)
