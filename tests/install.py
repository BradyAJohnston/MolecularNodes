import subprocess
import sys
import os


def main():
    python = os.path.realpath(sys.executable)

    commands = [
        f"{python} -m pip install -r requirements.txt",
        # test
        # f'{python} -m pip uninstall pytest-snapshot'
        f"{python} -m pip install pytest pytest-cov syrupy",
    ]

    for command in commands:
        subprocess.run(command.split(" "))


if __name__ == "__main__":
    main()
