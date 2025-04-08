import sys
import pytest

argv = sys.argv
argv = argv[argv.index("--") + 1 :]


# run this script like this:
# /Applications/Blender.app/Contents/MacOS/Blender -b -P tests/run.py -- . -v
# /Applications/Blender.app/Contents/MacOS/Blender -b -P tests/run.py -- . -k test_color_lookup_supplied


def main():
    # run the test suite, and we have to manually return the result value if non-zero
    # value is returned for a failing test
    if len(argv) == 0:
        result = pytest.main()
    else:
        result = pytest.main(argv)
    if result.value != 0:
        sys.exit(result.value)


if __name__ == "__main__":
    main()
