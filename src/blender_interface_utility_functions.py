import bpy

def create_new_collection(name, coll_parent = None):
    """
    Function for creating a new collection, that is linked to the given collection, 
    or the master collection if none is supplied.
    
    """

    coll_new = bpy.data.collections.get(name)

    if not coll_new:
        coll_new = bpy.data.collections.new(name)

        if coll_parent:
            coll_parent.children.link(coll_new)
        else:
            bpy.context.scene.collection.children.link(coll_new)
        
    bpy.context.view_layer.active_layer_collection = bpy.context.view_layer.layer_collection.children[name]

    return coll_new

col1 = create_new_collection("MolecularNodes")

col2 = create_new_collection("newcol", col1)
create_new_collection("newnewcol", col2)
