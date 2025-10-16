import site
import subprocess
import sys
import tomllib  # Python 3.11+ (Blender 4.0+)
from pathlib import Path


def install_dependency_group(pyproject_path, group_name):
    """
    Install packages from a specific dependency group in pyproject.toml.

    Args:
        pyproject_path: Path to pyproject.toml file
        group_name: Name of the dependency group to install
    """
    python_exe = sys.executable
    pyproject_file = Path(pyproject_path)

    if not pyproject_file.exists():
        print(f"Error: {pyproject_path} not found!")
        return

    print(f"Reading {pyproject_path}")
    print(f"Python executable: {python_exe}\n")

    # Read pyproject.toml
    with open(pyproject_file, "rb") as f:
        pyproject_data = tomllib.load(f)

    # Get the dependency group
    try:
        # Check for dependency-groups (PEP 735 format)
        if "dependency-groups" in pyproject_data:
            dependencies = pyproject_data["dependency-groups"].get(group_name, [])
        # Check for project.optional-dependencies
        elif (
            "project" in pyproject_data
            and "optional-dependencies" in pyproject_data["project"]
        ):
            dependencies = pyproject_data["project"]["optional-dependencies"].get(
                group_name, []
            )
        # Check for tool.poetry.group format
        elif "tool" in pyproject_data and "poetry" in pyproject_data["tool"]:
            poetry_groups = pyproject_data["tool"]["poetry"].get("group", {})
            dependencies = poetry_groups.get(group_name, {}).get("dependencies", {})
            # Convert poetry dict format to list
            if isinstance(dependencies, dict):
                dependencies = list(dependencies.keys())
        else:
            print("Error: Could not find dependency groups in pyproject.toml")
            return

    except KeyError:
        print(f"Error: Group '{group_name}' not found in pyproject.toml")
        return

    if not dependencies:
        print(f"No dependencies found in group '{group_name}'")
        return

    print(f"Found {len(dependencies)} packages in group '{group_name}':")
    for dep in dependencies:
        print(f"  - {dep}")
    print()

    # Install each package
    for package in dependencies:
        print(f"Installing {package}...")
        try:
            subprocess.check_call(
                [
                    python_exe,
                    "-m",
                    "pip",
                    "install",
                    package,
                    "--upgrade",
                    "--target",
                    site.getsitepackages()[0],
                ]
            )
            print(f"✓ {package} installed successfully\n")
        except subprocess.CalledProcessError as e:
            print(f"✗ Failed to install {package}: {e}\n")


# Configuration
PYPROJECT_PATH = "pyproject.toml"  # Change this to your file path
GROUP_NAME = "test"  # Change this to your group name (e.g., "dev", "test", "docs")

if __name__ == "__main__":
    install_dependency_group(PYPROJECT_PATH, GROUP_NAME)
