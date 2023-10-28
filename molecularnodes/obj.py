import bpy
import numpy as np

def create_object(name: str, collection: bpy.types.Collection, locations, bonds=[]) -> bpy.types.Object:
    """
    Create a mesh with the given name in the given collection, using the supplied
    vertex locations and, if provided, bonds as edges.

    Parameters
    ----------
    name : str
        The name of the mesh object to be created.
    collection : bpy.types.Collection
        The collection to which the mesh object will be added.
    locations : array-like
        The list of vertex locations for the mesh, an nx3 np array of locations.
        Each element in the list represents a 3D point (x, y, z) for a vertex.
    bonds : list of tuples, optional
        The list of vertex index pairs representing bonds as edges for the mesh.
        Each tuple should contain two vertex indices (e.g., (index1, index2)).

    Returns
    -------
    bpy.types.Object
        The newly created mesh object.

    Notes
    -----
    - The 'name' should be a unique identifier for the created mesh object.
    - The 'locations' list should contain at least three 3D points to define a 3D triangle.
    - If 'bonds' are not provided, the mesh will have no edges.
    - If 'bonds' are provided, they should be valid vertex indices within the 'locations' list.

    Example
    -------
    ```python
    # Create a mesh object named "MyMesh" in the collection "MyCollection"
    # with vertex locations and bond edges.
    locations = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]
    bonds = [(0, 1), (1, 2), (2, 0)]
    my_object = create_object("MyMesh", bpy.data.collections['Collection'], locations, bonds)
    ```
    """
    # create a new mesh
    MN_mesh = bpy.data.meshes.new(name)
    MN_mesh.from_pydata(locations, bonds, faces=[])
    MN_object = bpy.data.objects.new(name, MN_mesh)
    MN_object['type'] = 'molecule'
    collection.objects.link(MN_object)
    return MN_object


def add_attribute(object: bpy.types.Object, name: str, data, type="FLOAT", domain="POINT"):
    """
    Add an attribute to the given object's geometry on the given domain.

    Parameters
    ----------
        object : bpy.types.Object
            The object to which the attribute will be added.
        name : str
            The name of the attribute.
        data : array-like
            The data to be assigned to the attribute. For "FLOAT_VECTOR" attributes, it should be a 1D array
            representing the vector data.
        type : str, optional, default: "FLOAT"
            The data type of the attribute. Possible values are "FLOAT", "FLOAT_VECTOR", "INT", or "BOOLEAN".
        domain : str, optional, default: "POINT"
            The domain to which the attribute is added. Possible values are "POINT" or other domains supported
            by the object.

    Returns
    -------
        Any
            The newly created attribute, which can be further manipulated or used in the 3D environment.

    Notes
    -----
        - The function supports adding both scalar and vector attributes.
        - The "FLOAT_VECTOR" attribute requires the input data to be a 1D array, and it will be reshaped internally
          to represent vectors with 3 components (x, y, z).
    """

    if type == "FLOAT_VECTOR":
        att = object.data.attributes.new(name, type, domain)
        # currently vectors have to be added as a 1d array. may change in the future
        # but currently must be reshaped then added as a 'vector' but supplying a 1d array
        vec_1d = data.reshape(len(data) * 3).copy(order = 'c')
        att.data.foreach_set('vector', vec_1d)
    else:
        att = object.data.attributes.new(name, type, domain)
        att.data.foreach_set('value', data.copy(order = 'c'))
    
    return att

def get_attribute(obj: bpy.types.Object, att_name='position') -> np.array:
    """
    Retrieve an attribute from the object as a NumPy array.

    Parameters
    ----------
    obj : bpy.types.Object
        The Blender object from which the attribute will be retrieved.
    att_name : str, optional
        The name of the attribute to retrieve. Default is 'position'.
        
    Returns
    -------
    np.array
        The attribute data as a NumPy array.
    
    Notes
    -----
    - This function retrieves the specified attribute from the object and returns it as a NumPy array.
    - The function assumes that the attribute data type is one of ['INT', 'FLOAT', 'BOOLEAN', 'FLOAT_VECTOR'].

    Example
    -------
    ```python
    # Assuming 'my_object' is a Blender object with an attribute named 'my_attribute'
    attribute_data = get_attribute(my_object, 'my_attribute')
    ```
    """

    # Get the attribute from the object's mesh
    att = obj.to_mesh().attributes[att_name]

    # Map attribute values to a NumPy array based on the attribute data type
    if att.data_type in ['INT', 'FLOAT', 'BOOLEAN']:
        # Define the mapping of Blender data types to NumPy data types
        d_type = {'INT': int, 'FLOAT': float, 'BOOLEAN': bool}
        # Convert attribute values to a NumPy array with the appropriate data type
        att_array = np.array(list(map(lambda x: x.value, att.data.values())), dtype=d_type.get(att.data_type))
    elif att.data_type == "FLOAT_VECTOR":
        # Convert attribute vectors to a NumPy array
        att_array = np.array(list(map(lambda x: x.vector, att.data.values())))
    else:
        # Unsupported data type, return an empty NumPy array
        att_array = np.array([])

    return att_array
