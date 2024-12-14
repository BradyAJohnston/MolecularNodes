import subprocess
import sys
import os


def main():
    python = os.path.realpath(sys.executable)

    commands = [
        # include extra dependencies
        [
            python,
            "-m",
            "pip",
            "install",
            "uv",
        ],
        [
            python,
            "-m",
            "uv",
            "pip",
            "install",
            "-r",
            "pyproject.toml",
            "--all-extras",
        ],
        [
            python,
            "-m",
            "uv",
            "pip",
            "install",
            "-e",
            ".",
        ],
    ]

    for command in commands:
        subprocess.run(command)


if __name__ == "__main__":
    main()
