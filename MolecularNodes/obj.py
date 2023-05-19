import bpy
import numpy as np

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

def get_attribute(obj: bpy.types.Object, 
                  att_name = 'position') -> np.array:
    """Retrieve Attribute from Object as Numpy Array
    """
    att = obj.to_mesh().attributes[att_name]
    if att.data_type in ['INT', 'FLOAT', 'BOOLEAN']:
        d_type = {
            'INT': int, 
            'FLOAT': float, 
            'BOOLEAN': bool
        }
        att_array = np.array(list(map(
            lambda x: x.value, 
            att.data.values()
        )), dtype = d_type.get(att.data_type))
    elif att.data_type == "FLOAT_VECTOR":
        att_array = np.array(list(map(
            lambda x: x.vector, 
            att.data.values()
        )))
    
    return att_array