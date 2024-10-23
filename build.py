import glob
import os
import subprocess
import sys
from dataclasses import dataclass
from typing import List, Union
import requests
import json
import re
from pathlib import Path
import warnings
import shutil
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor

import bpy


def run_python(args: str | List[str]):
    python = os.path.realpath(sys.executable)

    if isinstance(args, str):
        args = [python] + args.split(" ")
    elif isinstance(args, list):
        args = [python] + args
    else:
        raise ValueError(
            "Arguments must be a string to split into individual arguments by space"
            "or a list of individual arguments already split"
        )

    subprocess.run(args)


try:
    import tomlkit  # type: ignore
except ModuleNotFoundError:
    run_python("-m pip install tomlkit")
    import tomlkit  # type: ignore

toml_path = "molecularnodes/blender_manifest.toml"
whl_path = "./molecularnodes/wheels"


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
    "biotite==0.41.2",
]


build_platforms = [
    windows_x64,
    linux_x64,
    macos_arm,
    macos_intel,
]


class Wheel:
    def __init__(self, url: str, log_path: str = None):
        self.url: str = url
        self.log_path = log_path
        self.platform = self.log_path_to_platform(self.log_path)

    @property
    def name(self):
        return self.url_to_name(self.url)

    @property
    def filepath(self):
        return f"molecularnodes/wheels/{self.name}"

    def file_exists(self):
        return os.path.exists(self.filepath)

    def url_to_name(self, url: str) -> str:
        return url.split("/")[-1]

    def log_path_to_platform(self, logpath: str) -> str:
        return logpath.removesuffix("_pip_log.txt")

    def download_from_url(self, skip_existing: bool = True):
        if skip_existing and self.file_exists():
            print(f"Skipping {self.name}, already downloaded.")
            return

        with requests.get(self.url, stream=True) as r:
            r.raise_for_status()
            with open(self.filepath, "wb") as f:
                print(f"Downloading {self.name}")
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)


class WheelHandler:
    def __init__(self):
        self.whls: List[Wheel] = []

    def add_urls_from_log(self, filepath: str):
        pattern = r"Downloading link ([^ ]+)"

        try:
            with open(filepath, "r") as file:
                for line in file:
                    match = re.search(pattern, line)
                    if match:
                        whl = Wheel(match.group(1), filepath)
                        print(f"{whl.filepath=}")
                        if whl.file_exists():
                            self.whls.append(whl)
        except Exception as e:
            print(f"Error reading {filepath}: {e}")

    def extract_log_files(self) -> None:
        for logpath in glob.glob("*_log.txt"):
            self.add_urls_from_log(logpath)

    def download_needed_wheels(self):
        def trigger_whl_download(whl: Wheel) -> None:
            whl.download_from_url()

        with ThreadPoolExecutor() as executor:
            result = executor.map(trigger_whl_download, self.whls)

            list(result)

    def write_to_json(self):
        dict_to_write = {"urls": [whl.url for whl in self.whls]}
        with open("whl_urls.json", "w") as file:
            json.dump(dict_to_write, file, indent=True)


def remove_whls():
    for whl_file in glob.glob(os.path.join(whl_path, "*.whl")):
        os.remove(whl_file)


def download_whls(
    platforms: Union[Platform, List[Platform]],
    required_packages: List[str] = required_packages,
    python_version="3.11",
    clean: bool = True,
):
    if isinstance(platforms, Platform):
        platforms = [platforms]

    if clean:
        remove_whls()

    for platform in platforms:
        run_python(
            f"-m pip download {' '.join(required_packages)} --dest ./molecularnodes/wheels --only-binary=:all: --python-version={python_version} --platform={platform.pypi_suffix} --log {platform.pypi_suffix}_pip_log.txt"
        )


def update_toml_whls(platforms):
    # Define the path for wheel files

    wheels_dir = "molecularnodes/wheels"
    wheel_files = glob.glob(f"{wheels_dir}/*.whl")
    wheel_files.sort()

    # Packages to remove
    packages_to_remove = {
        "pyarrow",
        "certifi",
        "charset_normalizer",
        "idna",
        "numpy",
        "requests",
        "urllib3",
    }

    # Filter out unwanted wheel files
    to_remove = []
    to_keep = []
    for whl in wheel_files:
        if any(pkg in whl for pkg in packages_to_remove):
            to_remove.append(whl)
        else:
            to_keep.append(whl)

    # Remove the unwanted wheel files from the filesystem
    for whl in to_remove:
        os.remove(whl)

    handler = WheelHandler()
    handler.extract_log_files()
    print(f"{handler.whls=}")
    handler.download_needed_wheels()
    handler.write_to_json()

    # Load the TOML file
    with open(toml_path, "r") as file:
        manifest = tomlkit.parse(file.read())

    # Update the wheels list with the remaining wheel files
    manifest["wheels"] = [f"./wheels/{os.path.basename(whl)}" for whl in to_keep]

    # Simplify platform handling
    if not isinstance(platforms, list):
        platforms = [platforms]
    manifest["platforms"] = [p.metadata for p in platforms]

    # Write the updated TOML file
    with open(toml_path, "w") as file:
        file.write(
            tomlkit.dumps(manifest)
            .replace('["', '[\n\t"')
            .replace("\\\\", "/")
            .replace('", "', '",\n\t"')
            .replace('"]', '",\n]')
        )


def clean_files(suffix: str = ".blend1") -> None:
    pattern_to_remove = f"molecularnodes/**/*{suffix}"
    for blend1_file in glob.glob(pattern_to_remove, recursive=True):
        os.remove(blend1_file)


def build_extension(split: bool = True) -> None:
    for suffix in [".blend1", ".MNSession"]:
        clean_files(suffix=suffix)

    if split:
        subprocess.run(
            f"{bpy.app.binary_path} --command extension build"
            " --split-platforms --source-dir molecularnodes --output-dir .".split(" ")
        )
    else:
        subprocess.run(
            f"{bpy.app.binary_path} --command extension build "
            "--source-dir molecularnodes --output-dir .".split(" ")
        )


def extract_download_urls(log_file_path):
    pattern = r"Downloading link ([^ ]+)"
    download_urls = []

    try:
        with open(log_file_path, "r") as file:
            for line in file:
                match = re.search(pattern, line)
                if match:
                    download_urls.append(match.group(1))
    except Exception as e:
        print(f"Error reading {log_file_path}: {e}")

    return download_urls


class PlatformLog:
    def __init__(self, path: Path):
        self.path = path
        self.urls = extract_download_urls(self.path)
        self.url_dict = self.urls_to_dict(self.urls)

    def urls_to_dict(self, urls: List[str]) -> dict:
        return {url.split("/")[-1]: url for url in urls}

    def remove_log(self) -> None:
        os.remove(self.path)
        del self


class LogHandler:
    def __init__(self, logs: List[str]):
        self.log_filenames = [Path(log) for log in logs]
        self.logs = [PlatformLog(x) for x in self.log_filenames]
        self.log_dict: dict = {}
        self.process_logs()

    def remove_logs(self):
        for log in self.logs:
            log.remove_log()


def build(platform) -> None:
    # download_whls(platform)
    update_toml_whls(platform)
    build_extension()


def logs_to_json():
    logs = glob.glob("*log.txt")
    if logs is None or len(logs) == 0:
        return None

    handler = LogHandler(logs)
    handler.process_logs()

    with open("download_urls.json", "w") as f:
        json.dump(handler.log_dict, f, indent=True)
    handler.remove_logs()


def main():
    # for platform in build_platforms:
    #     build(platform)
    build(build_platforms)
    # logs_to_json()


if __name__ == "__main__":
    main()
