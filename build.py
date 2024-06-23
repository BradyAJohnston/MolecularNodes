import glob
import os
import shutil
import subprocess
import sys
import zipfile
from dataclasses import dataclass
from typing import List, Union
import re

import tomlkit


@dataclass
class Platform:
    pypi_suffix: str
    metadata: str


# tags for blender metadata
# platforms = ["windows-x64", "macos-arm64", "linux-x64", "windows-arm64", "macos-x64"]


windows_x64 = Platform(pypi_suffix="win_amd64", metadata="windows-x64")
linux_x64 = Platform(pypi_suffix="manylinux2014_x86_64", metadata="linux-x64")
macos_arm = Platform(pypi_suffix="macosx_12_0_arm64", metadata="macos-arm64")
macos_intel = Platform(pypi_suffix="macosx_10_16_x86_64", metadata="macos-x64")


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
        for package in required_packages:
            run_python(
                f"-m pip download {package} --dest ./molecularnodes/wheels --only-binary=:all: --python-version={python_version} --platform={platform.pypi_suffix}"
            )


toml_path = "molecularnodes/blender_manifest.toml"


def get_version(string):
    str_match = r"\d+\.\d+.\d"
    return re.search(str_match, string).group()


def update_toml_whls(platform: Platform):
    # List all .whl files in the wheels/ subdirectory
    wheel_files = glob.glob("./wheels/*.whl", root_dir="molecularnodes")

    # scipy and pyarrow both are extremly large packages and aren't actually required by us,
    # so we can safely remove them (and potentially others) for bundling a smaller wheels/*
    for whl in wheel_files:
        packages_to_remove = [
            "pyarrow",
            "certifi",
            "charset_normalizer",
            "idna",
            "numpy",
            "requests",
            "urllib3",
        ]
        for name in packages_to_remove:
            if name in whl:
                whl_path = os.path.join("molecularnodes/", whl.removeprefix("./"))
                os.remove(whl_path)
                wheel_files.remove(whl)

    # Load the TOML file
    with open(toml_path, "r") as file:
        manifest = tomlkit.parse(file.read())

    # Update the wheels list
    manifest["wheels"] = wheel_files
    manifest["version"] = "{}-{}".format(
        get_version(manifest["version"]), platform.metadata
    )
    manifest["platforms"] = [platform.metadata]

    manifest_str = (
        tomlkit.dumps(manifest)
        .replace('["', '[\n\t"')
        .replace("\\\\", "/")
        .replace('", "', '",\n\t"')
        .replace('"]', '",\n]')
    )

    # Write the updated TOML file
    with open(toml_path, "w") as file:
        file.write(manifest_str)


def reset_toml() -> None:
    with open(toml_path, "r") as file:
        manifest = tomlkit.parse(file.read())

    manifest["version"] = get_version(manifest["version"])
    manifest["wheels"] = []
    manifest["platforms"] = []

    with open(toml_path, "w") as file:
        file.write(tomlkit.dumps(manifest))


def zip_extension(platform: Platform) -> None:
    # Load the TOML file
    with open(toml_path, "r") as file:
        manifest = tomlkit.parse(file.read())

    # Get the version number, which would have the platform appended onto it
    # for bundling for the blender extension platform and remove the platform to
    # keep just the version number

    zip_file_name = "molecularnodes-{}-{}.zip".format(
        get_version(manifest["version"]), platform.metadata
    )

    # Create the zip file
    with zipfile.ZipFile(zip_file_name, "w", zipfile.ZIP_DEFLATED) as zip_file:
        # Walk the molecularnodes folder
        for root, dirs, files in os.walk("molecularnodes"):
            # Ignore __pycache__ directories
            if "__pycache__" in root or "blend1" in files:
                continue
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
    reset_toml()


if __name__ == "__main__":
    for platform in build_platforms:
        build(platform)
