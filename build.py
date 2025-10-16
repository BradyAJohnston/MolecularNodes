import argparse
import glob
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from typing import List, Union
import tomlkit


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


TOML_PATH = "molecularnodes/blender_manifest.toml"
WHL_PATH = "./molecularnodes/wheels"
PYPROJ_PATH = "./pyproject.toml"
UVLOCK_PATH = "./uv.lock"


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


with open(PYPROJ_PATH, "r") as file:
    pyproj = tomlkit.parse(file.read())
    required_packages = pyproj["project"]["dependencies"]


build_platforms = [
    windows_x64,
    linux_x64,
    macos_arm,
    macos_intel,
]


def remove_whls():
    for whl_file in glob.glob(os.path.join(WHL_PATH, "*.whl")):
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
            f"-m pip download {' '.join(required_packages)} --dest ./molecularnodes/wheels --only-binary=:all: --python-version={python_version} --platform={platform.pypi_suffix}"
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
    with open(TOML_PATH, "r") as file:
        manifest = tomlkit.parse(file.read())

    # Update the wheels list with the remaining wheel files
    manifest["wheels"] = [f"./wheels/{os.path.basename(whl)}" for whl in to_keep]

    # Simplify platform handling
    if not isinstance(platforms, list):
        platforms = [platforms]
    manifest["platforms"] = [p.metadata for p in platforms]

    # Write the updated TOML file
    with open(TOML_PATH, "w") as file:
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


def build_extension(split: bool = True, blender_path: str = None) -> None:
    """Build the Blender extension.

    Args:
        split: Whether to build separate packages for each platform
        blender_path: Path to Blender executable. If None, tries to get from bpy module.
    """
    # Try to find Blender executable
    if not blender_path:
        try:
            import bpy

            blender_path = bpy.app.binary_path
        except (ImportError, AttributeError):
            pass

    if not blender_path:
        print("\nWarning: Blender executable path not available")
        print("The extension files have been prepared but not built.")
        print("To build the extension, either:")
        print(
            "  1. Run: blender --command extension build --split-platforms --source-dir molecularnodes --output-dir ."
        )
        print("  2. Re-run with: --blender-path /path/to/blender")
        return

    print(f"\nBuilding extension using Blender at: {blender_path}")

    for suffix in [".blend1", ".MNSession"]:
        clean_files(suffix=suffix)

    if split:
        subprocess.run(
            f"{blender_path} --command extension build"
            " --split-platforms --source-dir molecularnodes --output-dir .".split(" ")
        )
    else:
        subprocess.run(
            f"{blender_path} --command extension build "
            "--source-dir molecularnodes --output-dir .".split(" ")
        )


def get_all_dependencies_from_lock(package_name: str = "molecularnodes") -> set:
    """Get all transitive dependencies for a package from uv.lock.

    Args:
        package_name: The root package to start from

    Returns:
        A set of all package names (including transitive dependencies)
    """
    with open(UVLOCK_PATH, "r") as f:
        content = f.read()

    lock_data = tomlkit.parse(content)

    # Build a dependency graph
    dep_graph = {}
    for package in lock_data.get("package", []):
        name = package.get("name", "")
        deps = package.get("dependencies", [])
        dep_names = [d.get("name", "") for d in deps if isinstance(d, dict)]
        dep_graph[name] = dep_names

    # BFS to get all transitive dependencies
    all_deps = set()
    to_visit = [package_name]
    visited = set()

    while to_visit:
        current = to_visit.pop(0)
        if current in visited:
            continue
        visited.add(current)

        if current in dep_graph:
            for dep in dep_graph[current]:
                if dep not in visited:
                    all_deps.add(dep)
                    to_visit.append(dep)

    return all_deps


def parse_uv_lock_for_packages(package_names: set = None) -> dict:
    """Parse uv.lock and extract wheel info for specific packages.

    Args:
        package_names: Set of package names to extract. If None, extracts all packages.

    Returns:
        Dict mapping package name to {version, wheels: {filename: url}}
    """
    with open(UVLOCK_PATH, "r") as f:
        content = f.read()

    lock_data = tomlkit.parse(content)
    package_info_map = {}

    for package in lock_data.get("package", []):
        name = package.get("name", "")

        # Skip if we're filtering and this package isn't in the list
        if package_names is not None and name not in package_names:
            continue

        version = package.get("version", "")
        wheels = package.get("wheels", [])

        if wheels:
            normalized_name = name.lower().replace("_", "-")
            package_info_map[normalized_name] = {
                "version": version,
                "name": name,  # Keep original name
                "wheels": {},
            }

            for wheel in wheels:
                url = wheel.get("url", "")
                if url:
                    filename = url.split("/")[-1]
                    package_info_map[normalized_name]["wheels"][filename] = url

    return package_info_map


def download_wheels_from_lock(
    platforms: Union[Platform, List[Platform]],
    clean: bool = True,
    max_workers: int = 8,
    packages_to_exclude: set = None,
) -> None:
    """Download wheels from uv.lock for specified platforms.

    Args:
        platforms: Platform or list of platforms to download wheels for
        clean: Whether to remove existing wheel files before downloading
        max_workers: Maximum number of parallel download threads
        packages_to_exclude: Set of package names to exclude from download
    """
    import urllib.request
    from concurrent.futures import ThreadPoolExecutor, as_completed

    if isinstance(platforms, Platform):
        platforms = [platforms]

    if clean:
        remove_whls()

    # Ensure wheels directory exists
    os.makedirs(WHL_PATH, exist_ok=True)

    # Get all dependencies for molecularnodes
    print("Resolving dependencies from uv.lock...")
    all_deps = get_all_dependencies_from_lock("molecularnodes")
    print(f"Found {len(all_deps)} dependencies")

    # Filter out excluded packages
    if packages_to_exclude:
        all_deps = all_deps - packages_to_exclude
        print(f"After exclusions: {len(all_deps)} packages to download")

    # Parse uv.lock for these packages
    package_info = parse_uv_lock_for_packages(all_deps)

    # Collect all wheels to download for the specified platforms (using set to avoid duplicates)
    wheels_to_download = set()
    for pkg_name, pkg_data in package_info.items():
        for filename, url in pkg_data["wheels"].items():
            # Skip PyPy wheels - Blender doesn't support them
            if "-pp3" in filename or "pypy" in filename:
                continue

            # Universal wheels work for all platforms - add once
            if "py3-none-any" in filename or "py2.py3-none-any" in filename:
                wheels_to_download.add((filename, url))
                continue

            # Check if this wheel matches any of our target platforms
            # Uses Blender's platform matching logic from bl_extension_ops.py
            for platform in platforms:
                matched = False

                # For macOS, match any compatible macOS version with the same architecture
                if "macos" in platform.metadata:
                    if "arm64" in platform.metadata and (
                        "macosx" in filename and "arm64" in filename
                    ):
                        matched = True
                    elif "x64" in platform.metadata and (
                        "macosx" in filename and "x86_64" in filename
                    ):
                        matched = True
                    # Also match universal2 wheels for both macOS architectures
                    elif "macosx" in filename and "universal2" in filename:
                        matched = True
                # For Linux, match any manylinux wheel with the correct architecture
                # Blender accepts manylinux1, manylinux2010, manylinux2014, manylinux_2_XX, etc.
                elif "linux" in platform.metadata:
                    # Extract architecture from pypi_suffix (e.g., "manylinux2014_x86_64" -> "x86_64")
                    arch = platform.pypi_suffix.split("_", 1)[-1]  # Get everything after first underscore
                    if "manylinux" in filename and filename.endswith(("_" + arch, "_" + arch + ".whl")):
                        matched = True
                # For Windows, use exact suffix matching
                elif "windows" in platform.metadata and platform.pypi_suffix in filename:
                    matched = True

                if matched:
                    wheels_to_download.add((filename, url))
                    break

    # Convert set back to list for iteration
    wheels_to_download = list(wheels_to_download)

    print(f"Total wheels to download: {len(wheels_to_download)}")
    print(f"Using {max_workers} parallel download threads\n")

    def download_wheel(filename: str, url: str) -> tuple[str, bool, str]:
        """Download a single wheel file. Returns (filename, success, message)."""
        dest_path = os.path.join(WHL_PATH, filename)

        try:
            urllib.request.urlretrieve(url, dest_path)
            return (filename, True, "Downloaded successfully")
        except Exception as e:
            return (filename, False, str(e))

    # Download wheels in parallel
    success_count = 0
    failed_count = 0

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all download tasks
        future_to_wheel = {
            executor.submit(download_wheel, filename, url): (filename, url)
            for filename, url in wheels_to_download
        }

        # Process completed downloads
        for future in as_completed(future_to_wheel):
            filename, success, message = future.result()
            if success:
                print(f"✓ {filename}")
                success_count += 1
            else:
                print(f"✗ {filename}: {message}")
                failed_count += 1

    print(f"\nDownload complete: {success_count} succeeded, {failed_count} failed")


def download_wheels_from_manifest(clean: bool = True, max_workers: int = 8) -> None:
    """Download wheels listed in blender_manifest.toml using URLs from uv.lock.

    Args:
        clean: Whether to remove existing wheel files before downloading
        max_workers: Maximum number of parallel download threads
    """
    import urllib.request
    from concurrent.futures import ThreadPoolExecutor, as_completed

    if clean:
        remove_whls()

    # Ensure wheels directory exists
    os.makedirs(WHL_PATH, exist_ok=True)

    # Load the manifest to get the list of wheels
    with open(TOML_PATH, "r") as file:
        manifest = tomlkit.parse(file.read())

    wheels_in_manifest = manifest.get("wheels", [])

    # Parse uv.lock to get download URLs for all packages
    package_info_map = parse_uv_lock_for_packages()

    # Build exact match map from all packages
    exact_match_map = {}
    for pkg_data in package_info_map.values():
        exact_match_map.update(pkg_data["wheels"])

    print(f"Found {len(exact_match_map)} wheel URLs in uv.lock")
    print(f"Need to download {len(wheels_in_manifest)} wheels from manifest")
    print(f"Using {max_workers} parallel download threads\n")

    def parse_wheel_filename(filename: str) -> tuple[str, str, str]:
        """Parse wheel filename to extract package name, version, and platform tags.

        Returns (package_name, version, platform_tags)
        """
        # Pattern: package-version-pythontag-abitag-platformtag.whl
        match = re.match(r"^(.+?)-([\d\.]+.*?)-(py\d+|cp\d+).*\.whl$", filename)
        if match:
            pkg_name = match.group(1).lower().replace("_", "-")
            version = match.group(2)
            # Everything after version
            platform_tags = filename[
                len(match.group(1)) + 1 + len(match.group(2)) + 1 : -4
            ]
            return pkg_name, version, platform_tags
        return "", "", ""

    def find_matching_wheel(filename: str) -> tuple[str, str, str]:
        """Find matching wheel URL from uv.lock, handling version mismatches.

        Returns (url, actual_filename, message)
        """
        # Try exact match first
        if filename in exact_match_map:
            return (exact_match_map[filename], filename, "exact match")

        # Parse the wheel filename
        pkg_name, requested_version, platform_tags = parse_wheel_filename(filename)

        if not pkg_name:
            return ("", "", "Could not parse wheel filename")

        # Look up package in uv.lock
        if pkg_name not in package_info_map:
            return ("", "", f"Package '{pkg_name}' not found in uv.lock")

        lock_version = package_info_map[pkg_name]["version"]
        lock_wheels = package_info_map[pkg_name]["wheels"]

        # Try to find a wheel with matching platform tags
        for lock_filename, url in lock_wheels.items():
            _, _, lock_platform_tags = parse_wheel_filename(lock_filename)
            if lock_platform_tags == platform_tags:
                if requested_version != lock_version:
                    return (
                        url,
                        lock_filename,
                        f"version mismatch: requested {requested_version}, using {lock_version}",
                    )
                else:
                    return (url, lock_filename, "platform match")

        return ("", "", f"No matching platform found in uv.lock version {lock_version}")

    def download_wheel(wheel_path: str) -> tuple[str, str, bool, str]:
        """Download a single wheel file. Returns (requested_filename, actual_filename, success, message)."""
        requested_filename = os.path.basename(wheel_path)

        url, actual_filename, match_message = find_matching_wheel(requested_filename)

        if not url:
            return (requested_filename, "", False, match_message)

        dest_path = os.path.join(WHL_PATH, actual_filename)

        try:
            urllib.request.urlretrieve(url, dest_path)
            return (requested_filename, actual_filename, True, match_message)
        except Exception as e:
            return (requested_filename, actual_filename, False, str(e))

    # Download wheels in parallel
    success_count = 0
    failed_count = 0
    version_mismatch_count = 0

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all download tasks
        future_to_wheel = {
            executor.submit(download_wheel, wheel_path): wheel_path
            for wheel_path in wheels_in_manifest
        }

        # Process completed downloads
        for future in as_completed(future_to_wheel):
            requested_filename, actual_filename, success, message = future.result()
            if success:
                if "version mismatch" in message:
                    print(f"✓ {actual_filename} (substituted for {requested_filename})")
                    version_mismatch_count += 1
                else:
                    print(f"✓ {actual_filename}")
                success_count += 1
            else:
                print(f"✗ {requested_filename}: {message}")
                failed_count += 1

    print(f"\nDownload complete: {success_count} succeeded, {failed_count} failed")
    if version_mismatch_count > 0:
        print(f"  ({version_mismatch_count} version substitutions made)")


def verify_wheels_exist(
    platforms: Union[Platform, List[Platform]], packages_to_exclude: set = None
) -> tuple[bool, list, list]:
    """Verify that all required wheels exist in the wheels directory.

    Args:
        platforms: Platform or list of platforms to check for
        packages_to_exclude: Set of package names to exclude from verification

    Returns:
        A tuple of (all_exist, missing_packages, existing_files)
    """
    if isinstance(platforms, Platform):
        platforms = [platforms]

    # Get all dependencies for molecularnodes
    all_deps = get_all_dependencies_from_lock("molecularnodes")

    # Filter out excluded packages
    if packages_to_exclude:
        all_deps = all_deps - packages_to_exclude

    # Parse uv.lock for these packages
    package_info = parse_uv_lock_for_packages(all_deps)

    # Collect all expected wheels for the specified platforms (set already avoids duplicates)
    expected_wheels = set()
    for pkg_name, pkg_data in package_info.items():
        for filename, url in pkg_data["wheels"].items():
            # Skip PyPy wheels - Blender doesn't support them
            if "-pp3" in filename or "pypy" in filename:
                continue

            # Universal wheels work for all platforms - add once
            if "py3-none-any" in filename or "py2.py3-none-any" in filename:
                expected_wheels.add(filename)
                continue

            # Check if this wheel matches any of our target platforms
            # Uses Blender's platform matching logic from bl_extension_ops.py
            for platform in platforms:
                matched = False

                # For macOS, match any compatible macOS version with the same architecture
                if "macos" in platform.metadata:
                    if "arm64" in platform.metadata and (
                        "macosx" in filename and "arm64" in filename
                    ):
                        matched = True
                    elif "x64" in platform.metadata and (
                        "macosx" in filename and "x86_64" in filename
                    ):
                        matched = True
                    # Also match universal2 wheels for both macOS architectures
                    elif "macosx" in filename and "universal2" in filename:
                        matched = True
                # For Linux, match any manylinux wheel with the correct architecture
                # Blender accepts manylinux1, manylinux2010, manylinux2014, manylinux_2_XX, etc.
                elif "linux" in platform.metadata:
                    # Extract architecture from pypi_suffix (e.g., "manylinux2014_x86_64" -> "x86_64")
                    arch = platform.pypi_suffix.split("_", 1)[-1]  # Get everything after first underscore
                    if "manylinux" in filename and filename.endswith(("_" + arch, "_" + arch + ".whl")):
                        matched = True
                # For Windows, use exact suffix matching
                elif "windows" in platform.metadata and platform.pypi_suffix in filename:
                    matched = True

                if matched:
                    expected_wheels.add(filename)
                    break

    # Check which wheels actually exist
    existing_wheels = set()
    if os.path.exists(WHL_PATH):
        for whl_file in glob.glob(os.path.join(WHL_PATH, "*.whl")):
            existing_wheels.add(os.path.basename(whl_file))

    # Find missing wheels
    missing_wheels = expected_wheels - existing_wheels

    return (len(missing_wheels) == 0, sorted(missing_wheels), sorted(existing_wheels))


def build(
    platform,
    use_lock: bool = True,
    skip_download: bool = False,
    blender_path: str = None,
) -> None:
    """Build the extension.

    Args:
        platform: Platform or list of platforms to build for
        use_lock: If True, use uv.lock as source of truth. If False, use pip download.
        skip_download: If True, skip download and verify wheels exist before building.
        blender_path: Path to Blender executable for building extension.
    """
    # Packages that Blender already provides
    packages_to_exclude = {
        "pyarrow",
        "certifi",
        "charset-normalizer",
        "idna",
        "numpy",
        "requests",
        "urllib3",
    }

    if skip_download:
        print("Verifying all required packages exist...")
        all_exist, missing, existing = verify_wheels_exist(
            platform, packages_to_exclude
        )

        if not all_exist:
            print(f"\n✗ Missing {len(missing)} required wheel(s):")
            for whl in missing[:10]:  # Show first 10
                print(f"  - {whl}")
            if len(missing) > 10:
                print(f"  ... and {len(missing) - 10} more")
            print("\nRun without --build-only to download missing packages")
            return
        else:
            print(f"✓ All {len(existing)} required wheels are present")
    else:
        if use_lock:
            download_wheels_from_lock(platform, packages_to_exclude=packages_to_exclude)
        else:
            download_whls(platform)

    update_toml_whls(platform)
    build_extension(blender_path=blender_path)


def main():
    parser = argparse.ArgumentParser(
        description="Build Molecular Nodes Blender extension"
    )
    parser.add_argument(
        "--use-pip",
        action="store_true",
        help="Use pip download instead of uv.lock (legacy behavior)",
    )
    parser.add_argument(
        "--no-clean",
        action="store_true",
        help="Don't clean existing wheel files before downloading",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=8,
        help="Number of parallel download threads (default: 8)",
    )
    parser.add_argument(
        "--download-only",
        action="store_true",
        help="Only download wheels from uv.lock, don't update manifest or build extension",
    )
    parser.add_argument(
        "--build-only",
        action="store_true",
        help="Skip download, verify all required packages exist, then update manifest and build",
    )
    parser.add_argument(
        "--blender-path",
        type=str,
        help="Path to Blender executable (required for building extension)",
    )

    args = parser.parse_args()

    if args.download_only:
        print("Mode: Download wheels from uv.lock (download only)")
        packages_to_exclude = {
            "pyarrow",
            "certifi",
            "charset-normalizer",
            "idna",
            "numpy",
            "requests",
            "urllib3",
        }
        download_wheels_from_lock(
            build_platforms,
            clean=not args.no_clean,
            max_workers=args.workers,
            packages_to_exclude=packages_to_exclude,
        )
    else:
        use_lock = not args.use_pip
        mode_str = "pip download" if args.use_pip else "uv.lock"
        if args.build_only:
            print("Mode: Build extension only (verifying packages first)")
        else:
            print(f"Mode: Build extension using {mode_str}")
        build(
            build_platforms,
            use_lock=use_lock,
            skip_download=args.build_only,
            blender_path=args.blender_path,
        )


if __name__ == "__main__":
    main()
