import argparse
import re
import subprocess
import sys
import tomllib
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

TOML_PATH = Path("molecularnodes/blender_manifest.toml")
WHL_PATH = Path("molecularnodes/wheels")
UVLOCK_PATH = Path("uv.lock")
NODES_BLEND_PATH = Path("molecularnodes/assets/nodes.blend")

# Packages that Blender already provides, so their wheels must not be bundled.
# Canonical (hyphenated) package names; compare via normalize_package_name().
PACKAGES_TO_EXCLUDE = {
    "pyarrow",
    "certifi",
    "charset-normalizer",
    "idna",
    "numpy",
    "requests",
    "urllib3",
}


@dataclass
class Platform:
    pypi_suffix: str
    metadata: str


# tags for blender metadata
# platforms = ["windows-x64", "macos-arm64", "linux-x64", "windows-arm64", "macos-x64"]


build_platforms = [
    Platform(pypi_suffix="win_amd64", metadata="windows-x64"),
    Platform(pypi_suffix="manylinux2014_x86_64", metadata="linux-x64"),
    Platform(pypi_suffix="macosx_12_0_arm64", metadata="macos-arm64"),
    # unsure if these are the correct platform tags for windows_arm?
    # a build for this currently only ends up at 56 mb so I am guessing there are
    # some python packages that are not available for the plato
    # Platform(pypi_suffix="win_arm64", metadata="windows-arm64")
    # dropped intel mac support as of Blender 5.0
    # macos_intel = Platform(pypi_suffix="macosx_10_16_x86_64", metadata="macos-x64")
]


def normalize_package_name(name: str) -> str:
    """Normalize a package name for comparison (lowercase, hyphenated)."""
    return name.lower().replace("_", "-")


def package_name_from_wheel(filename: str) -> str:
    """Extract the normalized package name from a wheel filename."""
    return normalize_package_name(filename.split("-")[0])


def remove_whls() -> None:
    for whl_file in WHL_PATH.glob("*.whl"):
        whl_file.unlink()


def replace_toml_array(text: str, key: str, values: list[str]) -> str:
    """Replace or insert a top-level ``key = [...]`` array in TOML text.

    Only the array itself is touched; every other byte of the manifest
    (comments, ordering, other sections) is left as-is.
    """
    block = f"{key} = [\n" + "".join(f'\t"{value}",\n' for value in values) + "]"
    pattern = re.compile(rf"^{re.escape(key)}[ \t]*=[ \t]*\[[^\]]*\]", re.MULTILINE)
    if pattern.search(text):
        return pattern.sub(lambda _match: block, text, count=1)
    # Key not present: insert before the first table header, or append at the end.
    section = re.search(r"^\[", text, re.MULTILINE)
    insert_at = section.start() if section else len(text)
    return text[:insert_at] + block + "\n\n" + text[insert_at:]


def update_toml_whls(platforms: Platform | list[Platform]) -> None:
    """Point the manifest's wheels/platforms arrays at the downloaded wheels.

    Wheels for packages Blender already provides are deleted from disk and
    left out of the manifest.
    """
    if isinstance(platforms, Platform):
        platforms = [platforms]

    to_keep = []
    for whl in sorted(WHL_PATH.glob("*.whl")):
        if package_name_from_wheel(whl.name) in PACKAGES_TO_EXCLUDE:
            whl.unlink()
        else:
            to_keep.append(whl)

    text = TOML_PATH.read_text()
    text = replace_toml_array(text, "platforms", [p.metadata for p in platforms])
    text = replace_toml_array(
        text, "wheels", [f"./wheels/{whl.name}" for whl in to_keep]
    )
    TOML_PATH.write_text(text)


def clean_files(suffix: str = ".blend1") -> None:
    for file in Path("molecularnodes").rglob(f"*{suffix}"):
        file.unlink()


def build_extension(
    split: bool = True, blender_path: str | None = None, clean: bool = True
) -> None:
    """Build the Blender extension.

    Args:
        split: Whether to build separate packages for each platform
        blender_path: Path to the Blender executable to use when this script is
            not already running inside Blender. Defaults to "blender" on PATH.
        clean: Whether to remove stray .blend1/.MNSession files before building.
    """
    if not NODES_BLEND_PATH.exists():
        raise FileNotFoundError(
            f"{NODES_BLEND_PATH} is missing. It is a built asset library and is "
            "required in the packaged extension. Build it first with: "
            "uv run -m nodebpy.assets build"
        )

    if clean:
        for suffix in (".blend1", ".MNSession"):
            clean_files(suffix=suffix)

    try:
        import bpy

        executable = bpy.app.binary_path
        print(f"\nBuilding extension using current Blender instance: {executable}")
    except ImportError:
        executable = blender_path or "blender"
        print(f"\nBuilding extension using Blender at: {executable}")

    args = [executable, "--command", "extension", "build"]
    if split:
        args.append("--split-platforms")
    args += ["--source-dir", "molecularnodes", "--output-dir", "."]

    subprocess.run(args, check=True)
    print("Extension built successfully")


def get_all_dependencies_from_lock(package_name: str = "molecularnodes") -> set:
    """Get all transitive dependencies for a package from uv.lock.

    Args:
        package_name: The root package to start from

    Returns:
        A set of all package names (including transitive dependencies)
    """
    with open(UVLOCK_PATH, "rb") as f:
        lock_data = tomllib.load(f)

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


def parse_uv_lock_for_packages(package_names: set | None = None) -> dict:
    """Parse uv.lock and extract wheel info for specific packages.

    Args:
        package_names: Set of package names to extract. If None, extracts all packages.

    Returns:
        Dict mapping package name to {version, wheels: {filename: url}}
    """
    with open(UVLOCK_PATH, "rb") as f:
        lock_data = tomllib.load(f)
    package_info_map = {}

    for package in lock_data.get("package", []):
        name = package.get("name", "")

        # Skip if we're filtering and this package isn't in the list
        if package_names is not None and name not in package_names:
            continue

        version = package.get("version", "")
        wheels = package.get("wheels", [])

        if wheels:
            normalized_name = normalize_package_name(name)
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


def select_best_wheels(
    package_info: dict, platforms: list[Platform]
) -> dict[tuple[str, str], tuple[str, str, int]]:
    """Select the best wheel per package per platform.

    Skips PyPy wheels, prefers platform-specific wheels over universal ones,
    and uses Blender's platform matching logic from bl_extension_ops.py.

    Args:
        package_info: Output of parse_uv_lock_for_packages()
        platforms: Platforms to match wheels against

    Returns:
        Dict mapping (package_name, platform_metadata) to (filename, url, priority)
    """
    best_wheels: dict[tuple[str, str], tuple[str, str, int]] = {}

    def consider(
        pkg_name: str, platform: Platform, filename: str, url: str, priority: int
    ):
        key = (pkg_name, platform.metadata)
        # Only update if this wheel has higher priority
        if key not in best_wheels or best_wheels[key][2] < priority:
            best_wheels[key] = (filename, url, priority)

    for pkg_name, pkg_data in package_info.items():
        for filename, url in pkg_data["wheels"].items():
            # Skip PyPy wheels - Blender doesn't support them
            if "-pp3" in filename or "pypy" in filename:
                continue

            # Universal wheels work for all platforms, but prefer
            # platform-specific wheels (lower priority)
            if "py3-none-any" in filename or "py2.py3-none-any" in filename:
                for platform in platforms:
                    consider(pkg_name, platform, filename, url, priority=0)
                continue

            # Check if this wheel matches any of our target platforms
            for platform in platforms:
                matched = False
                priority = 1  # Default priority for platform-specific wheels

                # For macOS, match any compatible macOS version with the same architecture
                if "macos" in platform.metadata:
                    if "universal2" in filename and "macosx" in filename:
                        # Universal2 wheels work for both arm64 and x64 on macOS
                        matched = True
                        priority = 10  # High priority - works for both architectures
                    elif "arm64" in platform.metadata and (
                        "macosx" in filename and "arm64" in filename
                    ):
                        matched = True
                        priority = 5  # Architecture-specific wheel
                    elif "x64" in platform.metadata and (
                        "macosx" in filename and "x86_64" in filename
                    ):
                        matched = True
                        priority = 5  # Architecture-specific wheel
                # For Linux, match any manylinux wheel with the correct architecture
                # Blender accepts manylinux1, manylinux2010, manylinux2014, manylinux_2_XX, etc.
                elif "linux" in platform.metadata:
                    # Extract architecture from pypi_suffix
                    # (e.g., "manylinux2014_x86_64" -> "x86_64")
                    arch = platform.pypi_suffix.split("_", 1)[-1]
                    if "manylinux" in filename and ("_" + arch in filename):
                        matched = True
                        priority = 5
                # For Windows, use exact suffix matching
                elif (
                    "windows" in platform.metadata and platform.pypi_suffix in filename
                ):
                    matched = True
                    priority = 5

                if matched:
                    consider(pkg_name, platform, filename, url, priority)

    return best_wheels


def required_wheels_from_lock(
    platforms: Platform | list[Platform], packages_to_exclude: set | None = None
) -> dict[tuple[str, str], tuple[str, str, int]]:
    """Resolve the set of wheels required for the given platforms from uv.lock."""
    if isinstance(platforms, Platform):
        platforms = [platforms]

    # Get all dependencies for molecularnodes
    all_deps = get_all_dependencies_from_lock("molecularnodes")

    # Filter out excluded packages
    if packages_to_exclude:
        excluded = {normalize_package_name(pkg) for pkg in packages_to_exclude}
        all_deps = {
            dep for dep in all_deps if normalize_package_name(dep) not in excluded
        }

    # Parse uv.lock for these packages
    package_info = parse_uv_lock_for_packages(all_deps)

    return select_best_wheels(package_info, platforms)


def download_wheels_from_lock(
    platforms: Platform | list[Platform],
    clean: bool = True,
    max_workers: int = 8,
    packages_to_exclude: set | None = None,
) -> None:
    """Download wheels from uv.lock for specified platforms.

    Args:
        platforms: Platform or list of platforms to download wheels for
        clean: Whether to remove existing wheel files before downloading
        max_workers: Maximum number of parallel download threads
        packages_to_exclude: Set of package names to exclude from download
    """
    if clean:
        remove_whls()

    # Ensure wheels directory exists
    WHL_PATH.mkdir(parents=True, exist_ok=True)

    print("Resolving dependencies from uv.lock...")
    best_wheels = required_wheels_from_lock(platforms, packages_to_exclude)

    # Convert to set of (filename, url) tuples, removing duplicates and priorities
    wheels_to_download = list(
        set((filename, url) for filename, url, _ in best_wheels.values())
    )

    print(f"Total wheels to download: {len(wheels_to_download)}")
    print(f"Using {max_workers} parallel download threads\n")

    def download_wheel(filename: str, url: str) -> tuple[str, bool, str]:
        """Download a single wheel file. Returns (filename, success, message)."""
        try:
            urllib.request.urlretrieve(url, WHL_PATH / filename)
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
    if failed_count:
        raise RuntimeError(f"{failed_count} wheel download(s) failed")


def verify_wheels_exist(
    platforms: Platform | list[Platform], packages_to_exclude: set | None = None
) -> tuple[bool, list, list]:
    """Verify that all required wheels exist in the wheels directory.

    Args:
        platforms: Platform or list of platforms to check for
        packages_to_exclude: Set of package names to exclude from verification

    Returns:
        A tuple of (all_exist, missing_packages, existing_files)
    """
    best_wheels = required_wheels_from_lock(platforms, packages_to_exclude)

    # Convert to set of filenames, removing duplicates and priorities
    expected_wheels = set(filename for filename, _, _ in best_wheels.values())

    existing_wheels = {whl.name for whl in WHL_PATH.glob("*.whl")}

    missing_wheels = expected_wheels - existing_wheels

    return (len(missing_wheels) == 0, sorted(missing_wheels), sorted(existing_wheels))


def build(
    platforms: Platform | list[Platform],
    skip_download: bool = False,
    clean: bool = True,
    max_workers: int = 8,
) -> None:
    """Download (or verify) the wheels, update the manifest and build the extension.

    Args:
        platforms: Platform or list of platforms to build for
        skip_download: If True, skip download and verify wheels exist before building.
        clean: Whether to remove existing wheels before downloading and stray
            files before building.
        max_workers: Maximum number of parallel download threads.
    """
    if skip_download:
        print("Verifying all required packages exist...")
        all_exist, missing, existing = verify_wheels_exist(
            platforms, PACKAGES_TO_EXCLUDE
        )

        if not all_exist:
            print(f"\n✗ Missing {len(missing)} required wheel(s):")
            for whl in missing[:10]:  # Show first 10
                print(f"  - {whl}")
            if len(missing) > 10:
                print(f"  ... and {len(missing) - 10} more")
            print("\nRun without --build-only to download missing packages")
            sys.exit(1)
        print(f"✓ All {len(existing)} required wheels are present")
    else:
        download_wheels_from_lock(
            platforms,
            clean=clean,
            max_workers=max_workers,
            packages_to_exclude=PACKAGES_TO_EXCLUDE,
        )

    update_toml_whls(platforms)
    build_extension(clean=clean)


def parse_blender_args():
    """Parse arguments when run through Blender -P script or plain Python.

    Blender's sys.argv format: [blender_executable, -b, -P, script_name, -- script_args...]
    """
    # Find the -- separator that marks script arguments
    try:
        separator_index = sys.argv.index("--")
        script_args = sys.argv[separator_index + 1 :]
    except ValueError:
        # No -- separator. Inside Blender the remaining arguments are
        # Blender's own; under plain Python they belong to this script.
        try:
            import bpy  # noqa: F401

            script_args = []
        except ImportError:
            script_args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description="Build Molecular Nodes Blender extension"
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

    return parser.parse_args(script_args)


def main():
    args = parse_blender_args()

    if args.download_only:
        print("Mode: Download wheels from uv.lock (download only)")
        download_wheels_from_lock(
            build_platforms,
            clean=not args.no_clean,
            max_workers=args.workers,
            packages_to_exclude=PACKAGES_TO_EXCLUDE,
        )
    else:
        if args.build_only:
            print("Mode: Build extension only (verifying packages first)")
        else:
            print("Mode: Build extension using uv.lock")
        build(
            build_platforms,
            skip_download=args.build_only,
            clean=not args.no_clean,
            max_workers=args.workers,
        )


if __name__ == "__main__":
    main()
