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


def remove_whls():
    for whl_file in glob.glob(os.path.join(whl_path, "*.whl")):
        os.remove(whl_file)


def download_to_file(name, url):
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(f"molecularnodes/wheels/{name}", "wb") as f:
            print(f"Downloading {name}")
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return name


def download_whls(
    platforms: Union[Platform, List[Platform]],
    required_packages: List[str] = required_packages,
    python_version="3.11",
    use_logged_urls: bool = True,
    clean: bool = False,
):
    if isinstance(platforms, Platform):
        platforms = [platforms]

    if clean:
        remove_whls()

    if use_logged_urls:
        url_file = Path("download_urls.json")
        try:
            with open(url_file, "r") as f:
                url_dict = json.loads(f.read())
                if len(url_dict) == 0:
                    raise ValueError(
                        "download_urls.json is empty, re-downloading .whl files"
                    )

                urls = {}
                for platform in platforms:
                    urls.update(url_dict[platform.pypi_suffix])

                with ThreadPoolExecutor() as executor:
                    results = executor.map(download_to_file, urls.keys(), urls.values())

                    # print(list(results))

            return None
        except (FileNotFoundError, ValueError):
            warnings.warn("download_urls.json not found, re-downloading .whl files")

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

    @property
    def platform(self) -> str:
        return self.path.name.removesuffix("_pip_log.txt")

    def remove_log(self) -> None:
        os.remove(self.path)
        del self


class LogHandler:
    def __init__(self, logs: List[str]):
        self.log_filenames = [Path(log) for log in logs]
        self.logs = [PlatformLog(x) for x in self.log_filenames]
        self.log_dict: dict = {}
        self.process_logs()

    def process_logs(self):
        for log in self.logs:
            self.log_dict[log.platform] = log.url_dict

    def remove_logs(self):
        for log in self.logs:
            log.remove_log()


def build(platform) -> None:
    download_whls(platform)
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
    logs_to_json()


if __name__ == "__main__":
    main()
