import bpy

def coll_mn():
    coll = bpy.data.collections.get('MolecularNodes')
    if not coll:
        coll = bpy.data.collections.new('MolecularNodes')
        bpy.context.scene.collection.children.link(coll)
    return coll


# check if a particular property already exists or not
def property_exists(prop_path, glob, loc):
    try:
        eval(prop_path, glob, loc)
        return True
    except:
        return False
    
