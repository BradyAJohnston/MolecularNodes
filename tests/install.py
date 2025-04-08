import os
import subprocess
import sys


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
            "--extra=test",
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
