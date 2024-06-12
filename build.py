import glob
import os
import shutil
import subprocess
import sys
import zipfile
from dataclasses import dataclass
from typing import List, Union

import tomlkit


@dataclass
class Platform:
    meatadata: str
    pypi_suffix: str


# tags for blender metadata
# platforms = ["windows-x64", "macos-arm64", "linux-x64"]
# Other supported platforms: "windows-arm64", "macos-x64"
windows_x64 = Platform("windows-x64", "win_amd64")
linux_x64 = Platform("manylinux2014_28_x86_64", "linux-x64")
macos_arm = Platform("macos-arm64", "macosx_12_0_arm64")
macos_intel = Platform("macosx_10_9_x86_64", "macos-x64")


# download the required .whl files for all platforms

required_packages = [
    "mrcfile==1.4.3",
    "starfile==0.5.6",
    "MDAnalysis==2.7.0",
    "biotite==0.40.0",
]


build_platforms = [
    windows_x64,
    linux_x64,
    macos_arm,
    macos_intel,
]


def run_python(args: str):
    python = os.path.realpath(sys.executable)
    subprocess.run([python] + args.split(" "))


def download_whls(
    python_version="3.11",
    packages: List[str] = required_packages,
    platforms: Union[Platform, List[Platform]] = build_platforms,
):
    if isinstance(platforms, Platform):
        platforms = [platforms]

    for platform in platforms:
        run_python(
            f"-m pip download {' '.join(packages)} --dest ./molecularnodes/wheels --only-binary=:all: --python-version={python_version} --platform={platform.pypi_suffix}"
        )


toml_path = "molecularnodes/blender_manifest.toml"


def update_toml_whls(platform: Platform | None = None):
    # List all .whl files in the wheels/ subdirectory
    wheel_files = glob.glob("./wheels/*.whl", root_dir="molecularnodes")

    # Load the TOML file
    with open(toml_path, "r") as file:
        manifest = tomlkit.parse(file.read())

    # Update the wheels list
    manifest["wheels"] = wheel_files
    manifest["version"] = "{}-{}".format(manifest["version"], platform.meatadata)
    manifest["platform"] = platform.meatadata

    manifest_str = (
        tomlkit.dumps(manifest)
        .replace('["', '[\n\t"')
        .replace('", "', '",\n\t"')
        .replace('"]', '",\n]')
    )

    # Write the updated TOML file
    with open(toml_path, "w") as file:
        file.write(manifest_str)


def reset_version() -> None:
    with open(toml_path, "r") as file:
        manifest = tomlkit.parse(file.read())

    manifest["version"] = manifest["version"].split("-")[0]

    with open(toml_path, "w") as file:
        file.write(tomlkit.dumps(manifest))


def zip_extension(platform: Platform | None = None) -> None:
    # Load the TOML file
    with open(toml_path, "r") as file:
        manifest = tomlkit.parse(file.read())

    # Get the version number
    version = manifest["version"].split("-")[0]

    if platform:
        # Define the zip file name
        zip_file_name = f"molecularnodes_{version}_{platform.meatadata}.zip"
    else:
        zip_file_name = f"molecularnodes_{version}.zip"

    # Create the zip file
    with zipfile.ZipFile(zip_file_name, "w", zipfile.ZIP_DEFLATED) as zip_file:
        # Walk the molecularnodes folder
        for root, dirs, files in os.walk("molecularnodes"):
            for file in files:
                # Get the file path
                file_path = os.path.join(root, file)

                # Add the file to the zip file
                zip_file.write(
                    file_path, arcname=os.path.relpath(file_path, "molecularnodes")
                )


def remove_whls():
    dir_path = "./molecularnodes/wheels"
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


def build(platform: Platform) -> None:
    download_whls(platforms=platform)
    update_toml_whls(platform=platform)
    zip_extension(platform=platform)
    remove_whls()
    reset_version()


if __name__ == "__main__":
    for platform in build_platforms:
        build(platform)
