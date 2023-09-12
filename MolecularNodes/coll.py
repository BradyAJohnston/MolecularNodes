import bpy

def mn() -> bpy.types.Collection:
    """
    Return the `MolecularNodes` Collection
    
    Returns
    -------
    coll : bpy.types.Collection
        The 'MolecularNodes' collection inside the Blender scene. If it doesn't 
        exist, it will be created.
    """
    coll = bpy.data.collections.get('MolecularNodes')
    if not coll:
        coll = bpy.data.collections.new('MolecularNodes')
        bpy.context.scene.collection.children.link(coll)
    return coll

def data(name="data") -> bpy.types.Collection:
    """
    A collection for storing MN related data objects.

    Parameters
    ----------
    name : str, optional
        Name of the data collection. Default is "data".

    Returns
    -------
    coll : bpy.types.Collection
        The collection for storing MN related data objects.
    """
    name = f"MN_{name}"
    coll = bpy.data.collections.get(name)
    if not coll:
        coll = bpy.data.collections.new(name)
        mn().children.link(coll)
        
        # Disable the view of the data collection
        bpy.context.view_layer.layer_collection.children['MolecularNodes'].children[name].exclude = True
    return coll


def frames(name="", parent=None, prefix="frames") -> bpy.types.Collection:
    """
    Create a Collection for Frames of a Trajectory

    Parameters
    ----------
    name : str, optional
        Name of the collection for the frames. Default is "".
    parent : bpy.types.Collection, optional
        A blender collection which will become the parent collection. 
        Default is the MolecularNodes collection if None.

    Returns
    -------
    coll_frames : bpy.types.Collection
        The newly created collection for frames.
    """
    coll_frames = bpy.data.collections.new(f"{prefix}_{name}")
    if not parent:
        mn().children.link(coll_frames)
    else:
        parent.children.link(coll_frames)
    
    return coll_frames


