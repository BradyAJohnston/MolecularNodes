import bpy
import numpy as np

def create_object(
    locations: np.ndarray, 
    edges = [], 
    faces = [], 
    name: str = 'NewObject', 
    collection: bpy.types.Collection = None
    ) -> bpy.types.Object:
    """
    Create a mesh with the given name in the given collection, using the supplied
    vertex locations and, if provided, edges.

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
    faces : list of tupls, options
        List of vertex groups, each group will become a distrinct face.

    Returns
    -------
    bpy.types.Object
        The newly created mesh object.

    Notes
    -----
    - The 'name' should be a unique identifier for the created mesh object.
    - The 'locations' list should contain at least three 3D points to define a 3D triangle.
    - If 'edges' are not provided, the mesh will have no edges.
    - If 'edges' are provided, they should be valid vertex indices within the 'locations' list.

    Example
    -------
    ```python
    # Create a mesh object named "MyMesh" in the collection "MyCollection"
    # with vertex locations and bond edges.
    locations = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]
    edges = [(0, 1), (1, 2), (2, 0)]
    my_object = create_object(lcoations, edges, name='MySimpleObject')
    ```
    """
    # Have to first create a mesh, then fill th emesh with data from verts, edges and faces
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(locations, edges, faces=[])
    
    # create a new object that has a name and links to the created mesh
    object = bpy.data.objects.new(name, mesh)
    if not collection:
        collection = bpy.data.collections['Collection']
    collection.objects.link(object)
    
    # assign various simple attributes / properties to the object
    object['type'] = 'molecule'
    return object


def add_attribute(object: bpy.types.Object, name: str, data: np.ndarray, type="FLOAT", domain="POINT", overwrite: bool = False):
    """
    Add an attribute to the given object's geometry on the given domain.

    Parameters
    ----------
        object : bpy.types.Object
            The object to which the attribute will be added.
        name : str
            The name of the attribute.
        data : array-like
            The data to be assigned to the attribute. "FLOAT_VECTOR" and "FLOAT_COLOR" entries should be of length 3 and 4 respectively. 
        type : str, optional, default: "FLOAT"
            The data type of the attribute. Possible values are "FLOAT", "FLOAT_VECTOR", "FLOAT_COLOR", "INT", or "BOOLEAN".
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
    att = object.data.attributes.get(name)
    if not att or not overwrite:
        att = object.data.attributes.new(name, type, domain)
    
    types = {
        'FLOAT_VECTOR': 'vector',
        'FLOAT_COLOR': 'color', 
        'QUATERNION': 'value',
        'INT': 'value', 
        'FLOAT': 'value',
        'BOOLEAN': 'value' 
        
    }
    
    att.data.foreach_set(types[type], data.reshape(-1).copy(order='c'))
    
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
    elif att.data_type == "FLOAT_COLOR":
        att_array = np.array(list(map(lambda x: x.color, att.data.values())))
    elif att.data_type == "QUATERNION":
        att_array = np.array(list(map(lambda x: x.value, att.data.values())))
    else:
        # Unsupported data type, return an empty NumPy array
        att_array = np.array([])

    return att_array


def set_position(object, locations: np.ndarray):
    """
    Update the vertex positions of a Blender object.

    Parameters
    ----------
    object : bpy.types.Object
        The Blender object whose vertex positions need to be updated.
    locations : numpy.ndarray, optional
        An array containing the new vertex positions. Default is an empty array.

    Returns
    -------
    None

    Raises
    ------
    TypeError
        If `object` is not of type `bpy.types.Object`.
        If `locations` is not of type `numpy.ndarray`.
    ValueError
        If the shape of `locations` is not (n, 3), where n is the number of vertices.
    AttributeError
        If the object's data block does not have a 'position' attribute.

    Notes
    -----
    The `locations` array should be of shape (n, 3), where n is the number of vertices.
    The `object` should have a data block containing a 'position' attribute.

    Example
    -------
    set_position(obj, np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
    """
    # Check if the input object is valid
    if not isinstance(object, bpy.types.Object):
        raise TypeError("Expected 'object' to be a bpy.types.Object")

    # Check if the input locations array is valid
    if not isinstance(locations, np.ndarray):
        raise TypeError("Expected 'locations' to be a numpy.ndarray")

    if locations.shape[1] != 3:
        raise ValueError("The 'locations' array should be of shape (n, 3)")

    # Check if the object has a 'position' attribute
    if 'position' not in object.data.attributes:
        raise AttributeError("The object's data block must have a 'position' attribute")

    pos = object.data.attributes['position']

    # Ensure the locations array is flattened and compatible with the 'vector' attribute
    pos.data.foreach_set('vector', locations.reshape(-1))

    # Update the object's data
    object.data.update()
