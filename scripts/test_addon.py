import sys
import os
import subprocess

# install package for running the testing framework
# currently using my own custom branch to install correct version for Intel MacOS otherwise
# it incorrectly grabs the ARM build and the test fails

python_exe = os.path.realpath(sys.executable)
subprocess.run(
    [python_exe, '-m', 'pip', 'install', 'git+https://github.com/bradyajohnston/blender-addon-tester.git']
    # [python_exe, '-m', 'pip', 'install', 'C:\\Users\\BradyJohnston\\Documents\\GitHub\\blender-git\\blender-addon-tester']
)

try:
    import blender_addon_tester as bat
except ModuleNotFoundError as e:
    print(e)
    sys.exit(1)

def main():    
    if len(sys.argv) > 1:
        addon = sys.argv[1]
    else:
        addon = "MolecularNodes"
    if len(sys.argv) > 2:
        blender_rev = sys.argv[2]
    else:
        blender_rev = "3.5"
    
    try:
        exit_val = bat.test_blender_addon(addon_path=addon, blender_revision=blender_rev)
    except Exception as e:
        print(e)
        exit_val = 1
    sys.exit(exit_val)

if __name__ == "__main__":
    main()