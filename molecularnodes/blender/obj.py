import bpy
import numpy as np

from . import coll
from . import nodes
from dataclasses import dataclass


@dataclass
class AttributeTypeInfo:
    dname: str
    dtype: type
    width: int


TYPES = {key: AttributeTypeInfo(*values) for key, values in {
    'FLOAT_VECTOR': ('vector', float, 3),
    'FLOAT_COLOR': ('color', float, 4),
    'QUATERNION': ('value', float, 4),
    'INT': ('value', int, 1),
    'FLOAT': ('value', float, 1),
    'BOOLEAN': ('value', bool, 1)
}.items()}


class SimpleObject:
    def __init__(self, locations=None, edges=None, faces=None, name='NewObject', collection=None):
        self.object = create_object(locations, edges, faces, name, collection)

    def set_attribute(
        self,
        data: np.ndarray,
        name='NewAttribute',
        type=None,
        domain='POINT',
        overwrite=True
    ):

        set_attribute(
            object=self.object,
            name=name,
            data=data,
            type=type,
            domain=domain,
            overwrite=overwrite
        )


class AttributeMismatchError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


def create_object(
    vertices: np.ndarray = [],
    edges: np.ndarray = [],
    faces: np.ndarray = [],
    name: str = 'NewObject',
    collection: bpy.types.Collection = None
) -> bpy.types.Object:
    """
    Create a new Blender object, initialised with locations for each vertex. 

    If edges and faces are supplied then these are also created on the mesh.

    Parameters
    ----------
        vertices : np.ndarray, optional
            The vertices of the vertices as a numpy array. Defaults to None.
        edges : np.ndarray, optional
            The edges of the object as a numpy array. Defaults to None.
        faces : np.ndarray, optional
            The faces of the object as a numpy array. Defaults to None.
        name : str, optional
            The name of the object. Defaults to 'NewObject'.
        collection : bpy.types.Collection, optional
            The collection to link the object to. Defaults to None.

    Returns
    -------
        bpy.types.Object
            The created object.
    """
    mesh = bpy.data.meshes.new(name)

    mesh.from_pydata(vertices=vertices, edges=edges, faces=faces)

    object = bpy.data.objects.new(name, mesh)

    if not collection:
        # Add the object to the scene if no collection is specified
        collection = bpy.data.collections['Collection']

    collection.objects.link(object)

    return object


def set_attribute(
    object: bpy.types.Object,
    name: str,
    data: np.ndarray,
    type=None,
    domain="POINT",
    overwrite: bool = False
) -> bpy.types.Attribute:
    """
    Adds and sets the values of an attribute on the object.

    Parameters
    ----------
    object : bpy.types.Object
        The Blender object.
    name : str
        The name of the attribute.
    data : np.ndarray
        The attribute data as a numpy array.
    type : str, optional
        The data type of the attribute. Defaults to "FLOAT". Possbible values are (
            'FLOAT_VECTOR', 'FLOAT_COLOR", 'QUATERNION', 'FLOAT', 'INT', 'BOOLEAN'
        )
    domain : str, optional
        The domain of the attribute. Defaults to "POINT". Currenlty only ('POINT', 'EDGE', 
        'FACE') have been tested.
    overwrite : bool, optional
        Whether to overwrite an existing attribute with the same name. Defaults to False.

    Returns
    -------
    bpy.types.Attribute
        The added attribute.
    """

    dtype = data.dtype
    shape = data.shape
    # if the datatype isn't specified, try to guess the datatype based on the
    # datatype of the ndarray. This should work but ultimately won't guess between
    # the quaternion and color datatype, so will just default to color
    if not type:
        if len(shape == 1):
            if isinstance(dtype, int):
                type = "INT"
            elif isinstance(dtype, float):
                type = "FLOAT"
            elif isinstance(dtype, bool):
                type = "BOOL"
        else:
            if shape[1] == 3:
                type = "FLOAT_VECTOR"
            elif shape[1] == 4:
                type == "FLOAT_COLOR"

    attribute = object.data.attributes.get(name)
    if not attribute or not overwrite:
        attribute = object.data.attributes.new(name, type, domain)

    if len(data) != len(attribute.data):
        raise AttributeMismatchError(
            f"Length of input data {len(data)=} is not equal to the size of the domain {len(attribute.data)=}"
        )

    # the 'foreach_set' requires a 1D array, regardless of the shape of the attribute
    # it also requires the order to be 'c' or blender might crash!!
    attribute.data.foreach_set(
        TYPES[type].dname, data.reshape(-1).copy(order='c'))

    return attribute


def get_attribute(object: bpy.types.Object, name='position') -> np.ndarray:
    """
    Get the attribute data from the object.

    Parameters:
        object (bpy.types.Object): The Blender object.
        name (str, optional): The name of the attribute. Defaults to 'position'.

    Returns:
        np.ndarray: The attribute data as a numpy array.
    """

    # Get the attribute and some metadata about it from the object
    att = object.data.attributes[name]
    n_att = len(att.data)
    data_type = TYPES[att.data_type]
    width = data_type.width

    # data to and from attributes has to be given and taken as a 1D array
    # we have the initialise the array first with the appropriate length, then we can
    # fill it with the given data using the 'foreach_get' method which is super fast C++
    # internal method
    arr = np.zeros(n_att * width, dtype=data_type.dtype)
    # it is currently not really consistent, but to get the values you need to use one of
    # the 'value', 'vector', 'color' etc from the types dict. This I could only figure
    # out through trial and error. I assume this might be changed / improved in the future
    att.data.foreach_get(data_type.dname, arr)

    # if the attribute should be 2D, reshape it before returning the numpy array
    if width > 1:
        return arr.reshape((n_att, width))
    else:
        return arr


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
        raise AttributeError(
            "The object's data block must have a 'position' attribute")

    pos = object.data.attributes['position']

    # Ensure the locations array is flattened and compatible with the 'vector' attribute
    pos.data.foreach_set('vector', locations.reshape(-1))

    # Update the object's data
    object.data.update()


def evaluate_using_mesh(object):
    """
    Evaluate the object using a debug object. Some objects can't currently have their
    Geometry Node trees evaluated (such as volumes), so we source the geometry they create
    into a mesh object, which can be evaluated and tested.

    Parameters
    ----------
    object : bpy.types.Object
        The object to be evaluated.

    Returns
    -------
    bpy.types.Object

    Notes
    -----
    Intended for debugging only.
    """
    # create mesh an object that contains a single vertex
    debug = create_object(vertices=np.zeros((1, 3), dtype=float))
    mod = nodes.get_mod(debug)
    mod.node_group = nodes.create_debug_group()
    mod.node_group.nodes['Object Info'].inputs['Object'].default_value = bpy.data.objects[object.name]

    # This is super important, otherwise the evaluated object will not be updated
    debug.update_tag()
    dg = bpy.context.evaluated_depsgraph_get()
    evaluated = debug.evaluated_get(dg)

    return evaluated


def create_data_object(transforms_array, collection=None, name='CellPackModel', world_scale=0.01, fallback=False):
    obj_data = bpy.data.objects.get(name)
    if obj_data and fallback:
        return obj_data

    # still requires a unique call TODO: figure out why
    transforms_array = np.unique(transforms_array)

    # TODO: this recalculating of chain_ids I don't like, need to figure out a better way
    # to handle this
    chain_ids = np.unique(transforms_array['chain_id'], return_inverse=True)[1]
    locations = transforms_array['translation'] * world_scale

    if not collection:
        collection = coll.data()

    obj_data = create_object(locations, collection=collection, name=name)
    set_attribute(obj_data, 'rotation',
                  transforms_array['rotation'], 'QUATERNION', 'POINT')
    set_attribute(obj_data, 'assembly_id',
                  transforms_array['assembly_id'], 'INT', 'POINT')
    set_attribute(obj_data, 'chain_id', chain_ids, 'INT', 'POINT')
    try:
        set_attribute(obj_data, 'transform_id',
                      transforms_array['transform_id'], 'INT', 'POINT')
    except ValueError:
        pass

    return obj_data
