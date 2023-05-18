import bpy

def create_object(name, collection, locations, bonds=[]):
    """
    Creates a mesh with the given name in the given collection, from the supplied
    values for the locations of vertices, and if supplied, bonds as edges.
    """
    # create a new mesh
    mol_mesh = bpy.data.meshes.new(name)
    mol_mesh.from_pydata(locations, bonds, faces=[])
    mol_object = bpy.data.objects.new(name, mol_mesh)
    collection.objects.link(mol_object)
    return mol_object

def add_attribute(object, name, data, type = "FLOAT", domain = "POINT", add = True):
    if not add:
        return None
    attribute = object.data.attributes.new(name, type, domain)
    attribute.data.foreach_set('value', data)