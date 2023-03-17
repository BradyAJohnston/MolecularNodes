import sys
try:
    import blender_addon_tester as BAT
except Exception as e:
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
        blender_rev = "3.4"
    
    try:
        exit_val = BAT.test_blender_addon(addon_path=addon, blender_revision=blender_rev)
    except Exception as e:
        print(e)
        exit_val = 1
    sys.exit(exit_val)

main()