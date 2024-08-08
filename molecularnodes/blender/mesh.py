import bpy
import numpy as np

from typing import Optional
from enum import Enum

from . import coll
from . import nodes
from dataclasses import dataclass


@dataclass
class AttributeTypeInfo:
    dname: str
    dtype: type
    width: int


TYPES = {
    key: AttributeTypeInfo(*values)
    for key, values in {
        "FLOAT_VECTOR": ("vector", float, [3]),
        "FLOAT_COLOR": ("color", float, [4]),
        "QUATERNION": ("value", float, [4]),
        "INT": ("value", int, [1]),
        "FLOAT": ("value", float, [1]),
        "INT32_2D": ("value", int, [2]),
        "FLOAT4X4": ("value", float, [4, 4]),
        "BOOLEAN": ("value", bool, [1]),
    }.items()
}


class AttributeMismatchError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


def centre(array: np.array):
    return np.mean(array, axis=0)


def centre_weighted(array: np.ndarray, weight: np.ndarray):
    return np.sum(array * weight.reshape((len(array), 1)), axis=0) / np.sum(weight)


class ObjectTracker:
    """
    A context manager for tracking new objects in Blender.

    This class provides a way to track new objects that are added to Blender's bpy.data.objects collection.
    It stores the current objects when entering the context and provides a method to find new objects that were added when exiting the context.

    Methods
    -------
    new_objects():
        Returns a list of new objects that were added to bpy.data.objects while in the context.
    """

    def __enter__(self):
        """
        Store the current objects and their names when entering the context.

        Returns
        -------
        self
            The instance of the class.
        """
        self.objects = list(bpy.context.scene.objects)
        return self

    def __exit__(self, type, value, traceback):
        pass

    def new_objects(self):
        """
        Find new objects that were added to bpy.data.objects while in the context.

        Use new_objects()[-1] to get the most recently added object.

        Returns
        -------
        list
            A list of new objects.
        """
        obj_names = list([o.name for o in self.objects])
        current_objects = bpy.context.scene.objects
        new_objects = []
        for obj in current_objects:
            if obj.name not in obj_names:
                new_objects.append(obj)
        return new_objects

    def latest(self):
        """
        Get the most recently added object.

        This method returns the most recently added object to bpy.data.objects while in the context.

        Returns
        -------
        bpy.types.Object
            The most recently added object.
        """
        return self.new_objects()[-1]


def create_object(
    vertices: np.ndarray = [],
    edges: np.ndarray = [],
    faces: np.ndarray = [],
    name: str = "NewObject",
    collection: bpy.types.Collection = None,
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
        collection = bpy.data.collections["Collection"]

    collection.objects.link(object)

    object["type"] = "molecule"

    return object


class AttributeDataType(Enum):
    FLOAT_VECTOR = "FLOAT_VECTOR"
    FLOAT_COLOR = "FLOAT_COLOR"
    QUATERNION = "QUATERNION"
    FLOAT = "FLOAT"
    INT = "INT"
    BOOLEAN = "BOOLEAN"
    FLOAT4X4 = "FLOAT4X4"


def store_named_attribute(
    obj: bpy.types.Object,
    name: str,
    data: np.ndarray,
    data_type: Optional[str] = None,
    domain: str = "POINT",
    overwrite: bool = True,
) -> bpy.types.Attribute:
    """
    Adds and sets the values of an attribute on the object.

    Parameters
    ----------
    obj : bpy.types.Object
        The Blender object.
    name : str
        The name of the attribute.
    data : np.ndarray
        The attribute data as a numpy array.
    type : str, optional
        The data type of the attribute. Defaults to None. Possible values are:
        'FLOAT_VECTOR', 'FLOAT_COLOR', 'FLOAT4X4', 'QUATERNION', 'FLOAT', 'INT', 'BOOLEAN'
    domain : str, optional
        The domain of the attribute. Defaults to 'POINT'. Currently, only 'POINT', 'EDGE',
        and 'FACE' have been tested.
    overwrite : bool, optional
        Whether to overwrite an existing attribute with the same name. Defaults to True.

    Returns
    -------
    bpy.types.Attribute
        The added attribute.
    """

    # if the datatype isn't specified, try to guess the datatype based on the
    # datatype of the ndarray. This should work but ultimately won't guess between
    # the quaternion and color datatype, so will just default to color
    if data_type is None:
        dtype = data.dtype
        shape = data.shape

        if len(shape) == 1:
            if np.issubdtype(dtype, np.int_):
                data_type = "INT"
            elif np.issubdtype(dtype, np.float_):
                data_type = "FLOAT"
            elif np.issubdtype(dtype, np.bool_):
                data_type = "BOOLEAN"
        elif len(shape) == 3 and shape[1:] == (4, 4):
            data_type = "FLOAT4X4"
        else:
            if shape[1] == 3:
                data_type = "FLOAT_VECTOR"
            elif shape[1] == 4:
                data_type = "FLOAT_COLOR"
            else:
                data_type = "FLOAT"
    # catch if the data_type still wasn't determined and report info about the data
    if data_type is None:
        data_type = "FLOAT"
        # raise ValueError(
        #     f"Unable to determine data type for {data}, {shape=}, {dtype=}"
        # )

    attribute = obj.data.attributes.get(name)  # type: ignore
    if not attribute or not overwrite:
        attribute = obj.data.attributes.new(name, data_type, domain)  # type: ignore

    if len(data) != len(attribute.data):
        raise AttributeMismatchError(
            f"Data length {len(data)}, dimensions {data.shape} does not equal the size of the target domain {domain}, len={len(attribute.data)=}"
        )

    # the 'foreach_set' requires a 1D array, regardless of the shape of the attribute
    # it also requires the order to be 'c' or blender might crash!!
    attribute.data.foreach_set(TYPES[data_type].dname, data.reshape(-1))

    # The updating of data doesn't work 100% of the time (see:
    # https://projects.blender.org/blender/blender/issues/118507) so this resetting of a
    # single vertex is the current fix. Not great as I can see it breaking when we are
    # missing a vertex - but for now we shouldn't be dealing with any situations where this
    # is the case For now we will set a single vert to it's own position, which triggers a
    # proper refresh of the object data.
    try:
        obj.data.vertices[0].co = obj.data.vertices[0].co  # type: ignore
    except AttributeError:
        obj.data.update()  # type: ignore

    return attribute


def named_attribute(
    object: bpy.types.Object, name="position", evaluate=False
) -> np.ndarray:
    """
    Get the attribute data from the object.

    Parameters:
        object (bpy.types.Object): The Blender object.
        name (str, optional): The name of the attribute. Defaults to 'position'.

    Returns:
        np.ndarray: The attribute data as a numpy array.
    """
    if evaluate:
        object = evaluated(object)
    attribute_names = object.data.attributes.keys()
    verbose = False
    if name not in attribute_names:
        if verbose:
            raise AttributeError(
                f"The selected attribute '{name}' does not exist on the mesh. \
                Possible attributes are: {attribute_names=}"
            )
        else:
            raise AttributeError(
                f"The selected attribute '{name}' does not exist on the mesh."
            )

    # Get the attribute and some metadata about it from the object
    att = object.data.attributes[name]
    n_att = len(att.data)
    data_type = TYPES[att.data_type]
    dim = data_type.width
    n_values = n_att
    for dimension in dim:
        n_values *= dimension

    # data to and from attributes has to be given and taken as a 1D array
    # we have the initialise the array first with the appropriate length, then we can
    # fill it with the given data using the 'foreach_get' method which is super fast C++
    # internal method
    array = np.zeros(n_values, dtype=data_type.dtype)
    # it is currently not really consistent, but to get the values you need to use one of
    # the 'value', 'vector', 'color' etc from the types dict. This I could only figure
    # out through trial and error. I assume this might be changed / improved in the future
    att.data.foreach_get(data_type.dname, array)

    if dim == [1]:
        return array
    else:
        # return an array with one row per item, even if a 1D attribute. Does this make sense?
        return array.reshape((n_att, *dim))


def import_vdb(file: str, collection: bpy.types.Collection = None) -> bpy.types.Object:
    """
    Imports a VDB file as a Blender volume object, in the MolecularNodes collection.

    Parameters
    ----------
    file : str
        Path to the VDB file.

    Returns
    -------
    bpy.types.Object
        A Blender object containing the imported volume data.
    """

    # import the volume object
    with ObjectTracker() as o:
        bpy.ops.object.volume_import(filepath=file, files=[])
        object = o.latest()

    if collection:
        # Move the object to the MolecularNodes collection
        initial_collection = object.users_collection[0]
        initial_collection.objects.unlink(object)
        collection = coll.mn()
        collection.objects.link(object)

    return object


def evaluated(object):
    "Return an object which has the modifiers evaluated."
    object.update_tag()
    return object.evaluated_get(bpy.context.evaluated_depsgraph_get())


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
    # create an empty mesh object. It's modifiers can be evaluated but some other
    # object types can't be currently through the API
    debug = create_object()
    mod = nodes.get_mod(debug)
    mod.node_group = nodes.create_debug_group()
    mod.node_group.nodes["Object Info"].inputs["Object"].default_value = object

    # need to use 'evaluate' otherwise the modifiers won't be taken into account
    return evaluated(debug)


def create_data_object(
    array, collection=None, name="DataObject", world_scale=0.01, fallback=False
):
    # still requires a unique call TODO: figure out why
    # I think this has to do with the bcif instancing extraction
    # array = np.unique(array)
    locations = array["translation"] * world_scale

    if not collection:
        collection = coll.data()

    object = create_object(locations, collection=collection, name=name)

    attributes = [
        ("rotation", "QUATERNION"),
        ("assembly_id", "INT"),
        ("chain_id", "INT"),
        ("transform_id", "INT"),
    ]

    for column, type in attributes:
        try:
            data = array[column]
        except ValueError:
            continue
        # us the unique sorted integer encoding version of the non-numeric
        # attribute, as GN doesn't support strings currently
        if np.issubdtype(data.dtype, str):
            data = np.unique(data, return_inverse=True)[1]

        store_named_attribute(
            object, name=column, data=data, data_type=type, domain="POINT"
        )

    return object
