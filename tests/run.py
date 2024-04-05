import pytest
import sys
argv = sys.argv
argv = argv[argv.index("--") + 1:]


# run this script like this:
# /Applications/Blender.app/Contents/MacOS/Blender -b -P tests/run.py -- . -v
# /Applications/Blender.app/Contents/MacOS/Blender -b -P tests/run.py -- . -k test_color_lookup_supplied

def main():
    pytest.main(argv)


if __name__ == "__main__":
    main()
