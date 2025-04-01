import os
import subprocess
import sys

argv = sys.argv
argv = argv[argv.index("--") + 1 :]


def main():
    python = os.path.realpath(sys.executable)
    subprocess.run([python] + argv)


if __name__ == "__main__":
    main()
