import subprocess
import sys
import os


def main():

    python = os.path.realpath(sys.executable)

    commands = [
        f'{python} -m pip install .',
        f'{python} -m pip install quartodoc'
    ]

    for command in commands:
        subprocess.run(command.split(' '))


if __name__ == "__main__":
    main()
