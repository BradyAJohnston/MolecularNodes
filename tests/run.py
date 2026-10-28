import sys
import pytest

argv = sys.argv
argv = argv[argv.index("--") + 1 :]


# run this script like this:
# /Applications/Blender.app/Contents/MacOS/Blender -b -P tests/run.py -- . -v
# /Applications/Blender.app/Contents/MacOS/Blender -b -P tests/run.py -- . -k test_color_lookup_supplied


def main():
    # bpy.ops.preferences.addon_enable(module="bl_ext.user_default.molecularnodes")
    # bpy.ops.wm.save_userpref()
    # run the test suite, and we have to manually return the result value if non-zero
    # value is returned for a failing test
    if len(argv) == 0:
        result = pytest.main()
    else:
        result = pytest.main(argv)
    # pytest.main returns an ExitCode, or a plain int from a plugin
    if int(result) != 0:
        sys.exit(int(result))


if __name__ == "__main__":
    main()
