import subprocess
import sys
import tempfile
from pathlib import Path


def install_locked(extra):
    """
    Install the versions pinned in uv.lock into Blender's Python, so tests in Blender
    run against the same dependencies as `uv sync`.

    Args:
        extra: Optional dependency group to install alongside the dependencies
    """
    python_exe = sys.executable
    subprocess.check_call([python_exe, "-m", "pip", "install", "uv"])
    uv = [python_exe, "-m", "uv"]
    with tempfile.TemporaryDirectory() as tmp:
        requirements = str(Path(tmp) / "requirements.txt")
        subprocess.check_call(
            [
                *uv,
                "export",
                "--quiet",
                "--frozen",
                "--extra",
                extra,
                "--no-emit-project",
                "--output-file",
                requirements,
            ]
        )
        subprocess.check_call(
            [*uv, "pip", "install", "--python", python_exe, "-r", requirements]
        )
    # the dependencies are already installed from the lockfile
    subprocess.check_call(
        [*uv, "pip", "install", "--python", python_exe, "--no-deps", "-e", "."]
    )


EXTRA = "dev"

if __name__ == "__main__":
    try:
        install_locked(EXTRA)
    except subprocess.CalledProcessError as e:
        print(f"Installation failed: {e}")
        sys.exit(e.returncode)
