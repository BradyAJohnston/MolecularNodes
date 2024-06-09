import subprocess
import sys
import os
import pathlib

REQUIREMENTS = pathlib.Path(
    pathlib.Path(__file__).resolve().parent.parent
) / "molecularnodes/requirements.txt"


def main():

    python = os.path.realpath(sys.executable)

    commands = [
        f'{python} -m pip install -r molecularnodes/requirements.txt',
        # f'{python} -m pip uninstall pytest-snapshot'
        f'{python} -m pip install pytest pytest-cov syrupy'
    ]

    for command in commands:
        subprocess.run(command.split(' '))


if __name__ == "__main__":
    main()
