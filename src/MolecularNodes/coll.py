import bpy

def mn():
    """Return the MolecularNodes Collection
    
    The collection called 'MolecularNodes' inside the Blender scene is returned. If the 
    collection does not exist first, it is created.
    """
    coll = bpy.data.collections.get('MolecularNodes')
    if not coll:
        coll = bpy.data.collections.new('MolecularNodes')
        bpy.context.scene.collection.children.link(coll)
    return coll

def frames(name="", parent=None, suffix="_frames"):
    """Create a Collection for Frames of a Trajectory

    Args:
        name (str, optional): Name of the collection for the frames. Defaults to "".
        parent (_type_, optional): A blender collection which will become the parent 
        collection. Defaults to the MolecularNodes collection if None.
    """
    coll_frames = bpy.data.collections.new(name + suffix)
    if not parent:
        mn().children.link(coll_frames)
    else:
        parent.children.link(coll_frames)
    
    return coll_frames

