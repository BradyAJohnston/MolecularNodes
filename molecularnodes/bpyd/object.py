import bpy
from bpy.types import Object
import numpy as np
from .attribute import (
    AttributeTypes,
    AttributeType,
    Domains,
    DomainType,
)
from . import attribute as attr
from .utils import centre
from mathutils import Matrix


class ObjectMissingError(Exception):
    def __init__(self, message: str):
        self.message = message
        super().__init__(self.message)


class ObjectTracker:
    """
    A context manager for tracking new objects in Blender.

    This class provides a way to track new objects that are added to Blender's bpy.data.objects collection.
    It stores the current objects when entering the context and provides a method to find new objects that were added when exiting the context.

    Methods
    -------
    new_objects()
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
        Object
            The most recently added object.
        """
        return self.new_objects()[-1]


class BlenderObject:
    """
    A convenience class for working with Blender objects
    """

    def __init__(self, obj: Object | None):
        """
        Initialize the BlenderObject.

        Parameters
        ----------
        obj : Object | None
            The Blender object to wrap.
        """
        if not isinstance(obj, Object):
            raise ValueError(f"{obj} must be a Blender object of type Object")
        self._object = obj

    @property
    def object(self) -> Object | None:
        """
        Get the Blender object.

        Returns
        -------
        Object | None
            The Blender object, or None if not found.
        """
        # If we don't have connection to an object, attempt to re-stablish to a new
        # object in the scene with the same UUID. This helps if duplicating / deleting
        # objects in the scene, but sometimes Blender just loses reference to the object
        # we are working with because we are manually setting the data on the mesh,
        # which can wreak havoc on the object database. To protect against this,
        # if we have a broken link we just attempt to find a new suitable object for it
        try:
            # if the connection is broken then trying to the name will raise a connection
            # error. If we are loading from a saved session then the object_ref will be
            # None and get an AttributeError
            self._object.name
            return self._object
        except (ReferenceError, AttributeError):
            for obj in bpy.data.objects:
                if obj.mn.uuid == self.uuid:
                    print(
                        Warning(
                            f"Lost connection to object: {self._object}, now connected to {obj}"
                        )
                    )
                    self._object = obj
                    return obj

            return None

    @object.setter
    def object(self, value: Object) -> None:
        """
        Set the Blender object.

        Parameters
        ----------
        value : Object
            The Blender object to set.
        """
        self._object = value

    def store_named_attribute(
        self,
        data: np.ndarray,
        name: str,
        atype: str | AttributeType | None = None,
        domain: str | DomainType = Domains.POINT,
    ) -> None:
        """
        Store a named attribute on the Blender object.

        Parameters
        ----------
        data : np.ndarray
            The data to be stored as an attribute.
        name : str
            The name for the attribute. Will overwrite an already existing attribute.
        atype : str or AttributeType or None, optional
            The attribute type to store the data as. Either string or selection from the
            AttributeTypes enum. None will attempt to infer the attribute type from the
            input array.
        domain : str or DomainType, optional
            The domain to store the attribute on. Defaults to Domains.POINT.

        Returns
        -------
        self
        """
        attr.store_named_attribute(
            self.object, data=data, name=name, atype=atype, domain=domain
        )
        return self

    def remove_named_attribute(self, name: str) -> None:
        """
        Remove a named attribute from the object.

        Parameters
        ----------
        name : str
            The name of the attribute to remove.
        """
        attr.remove_named_attribute(self.object, name=name)

    def named_attribute(self, name: str, evaluate: bool = False) -> np.ndarray:
        """
        Retrieve a named attribute from the object.

        Optionally, evaluate the object before reading the named attribute

        Parameters
        ----------
        name : str
            Name of the attribute to get.
        evaluate : bool, optional
            Whether to evaluate the object before reading the attribute (default is False).
        Returns
        -------
        np.ndarray
            The attribute read from the mesh as a numpy array.
        """
        return attr.named_attribute(self.object, name=name, evaluate=evaluate)

    def set_boolean(self, array: np.ndarray, name: str) -> None:
        """
        Store a boolean attribute on the Blender object.

        Parameters
        ----------
        array : np.ndarray
            The boolean data to be stored as an attribute.
        name : str
            The name for the attribute.
        """
        self.store_named_attribute(array, name=name, atype=AttributeTypes.BOOLEAN)

    def evaluate(self):
        """
        Evaluate the object and return a new BlenderObject with the evaluated object.

        Returns
        -------
        BlenderObject
            A new BlenderObject with the evaluated object.
        """
        obj = self.object
        obj.update_tag()
        evluated_obj = obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
        return BlenderObject(evluated_obj)

    def centroid(self, weight: str | np.ndarray | None = None) -> np.ndarray:
        """
        Return the centroid, potentially weighted by an attribute.

        If the weight is a string, an attribute of that name is attempted to be accessed
        on the mesh. If an array is given that array is used as weights. A value of None
        returns just the centroid calculation.

        Parameters
        ----------
        weight : str | np.ndarray | None, optional
            The weights to apply to the positions when calculating the centroid. Defaults to None.

        Returns
        -------
        np.ndarray
            A 3-component vector with the calculated centroid.
        """
        if isinstance(weight, str):
            return centre(self.position, self.named_attribute(weight))

        if isinstance(weight, np.ndarray):
            return centre(self.position, weight)

        return centre(self.position)

    @property
    def attributes(self):
        """
        Get the attributes of the Blender object.

        Returns
        -------
        bpy.types.Attributes
            The attributes of the Blender object.
        """
        return self.object.data.attributes

    @property
    def vertices(self):
        """
        Get the vertices of the Blender object.

        Returns
        -------
        bpy.types.Vertices
            The vertices of the Blender object.
        """
        return self.object.data.vertices

    @property
    def edges(self):
        """
        Get the edges of the Blender object.

        Returns
        -------
        bpy.types.Edges
            The edges of the Blender object.
        """
        return self.object.data.edges

    def transform_origin(self, matrix: Matrix) -> None:
        """
        Transform the origin of the Blender object.

        Parameters
        ----------
        matrix : Matrix
            The transformation matrix to apply to the origin.
        """
        self.object.matrix_local = matrix * self.object.matrix_world

    def transform_points(self, matrix: Matrix) -> None:
        """
        Transform the points of the Blender object.

        Parameters
        ----------
        matrix : Matrix
            The transformation matrix to apply to the points.
        """
        self.position = self.position * matrix

    @property
    def selected(self) -> np.ndarray:
        """
        Get the selected vertices of the Blender object.

        Returns
        -------
        np.ndarray
            The selected vertices of the Blender object.
        """
        return self.named_attribute(".select_vert")

    @property
    def name(self) -> str:
        """
        Get the name of the Blender object.

        Returns
        -------
        str
            The name of the Blender object.
        """
        obj = self.object
        if obj is None:
            return None

        return obj.name

    @name.setter
    def name(self, value: str) -> None:
        """
        Set the name of the Blender object.

        Parameters
        ----------
        value : str
            The name to set for the Blender object.
        """
        obj = self.object
        if obj is None:
            raise ObjectMissingError
        obj.name = value

    @property
    def position(self) -> np.ndarray:
        """
        Get the position of the vertices of the Blender object.

        Returns
        -------
        np.ndarray
            The position of the vertices of the Blender object.
        """
        return self.named_attribute("position")

    @position.setter
    def position(self, value: np.ndarray) -> None:
        """
        Set the position of the vertices of the Blender object.

        Parameters
        ----------
        value : np.ndarray
            The position to set for the vertices of the Blender object.
        """
        self.store_named_attribute(
            value,
            name="position",
            atype=AttributeTypes.FLOAT_VECTOR,
            domain=Domains.POINT,
        )

    def selected_positions(self, mask: np.ndarray | None = None) -> np.ndarray:
        """
        Get the positions of the selected vertices, optionally filtered by a mask.

        Parameters
        ----------
        mask : np.ndarray | None, optional
            The mask to filter the selected vertices. Defaults to None.

        Returns
        -------
        np.ndarray
            The positions of the selected vertices.
        """
        if mask is not None:
            return self.position[np.logical_and(self.selected, mask)]

        return self.position[self.selected]

    def list_attributes(
        self, evaluate: bool = False, drop_hidden: bool = False
    ) -> list | None:
        """
        Returns a list of attribute names for the object.

        Parameters
        ----------
        evaluate : bool, optional
            Whether to first evaluate the modifiers on the object before listing the
            available attributes.
        drop_hidden : bool, optional
            Whether to drop hidden attributes (those starting with a dot). Defaults to False.

        Returns
        -------
        list[str] | None
            A list of attribute names if the molecule object exists, None otherwise.
        """
        if evaluate:
            strings = list(self.evaluate().attributes.keys())
        else:
            strings = list(self.object.attributes.keys())

        if not drop_hidden:
            return strings
        else:
            return filter(lambda x: not x.startswith("."), strings)

    def __len__(self) -> int:
        """
        Get the number of vertices in the Blender object.

        Returns
        -------
        int
            The number of vertices in the Blender object.
        """
        return len(self.object.data.vertices)


def create_object(
    vertices: np.ndarray | None = None,
    edges: np.ndarray | None = None,
    faces: np.ndarray | None = None,
    name: str = "NewObject",
    collection: bpy.types.Collection | None = None,
) -> Object:
    """
    Create a new Blender object and corresponding mesh.

    Vertices are created for each row in the vertices array. If edges and / or faces are created then they are also
    initialized but default to None.

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
        Object
            The created object.
    """
    if vertices is None:
        vertices = []
    if edges is None:
        edges = []
    if faces is None:
        faces = []
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(vertices=vertices, edges=edges, faces=faces)
    obj = bpy.data.objects.new(name, mesh)
    if not collection:
        collection = bpy.data.collections["Collection"]
    collection.objects.link(obj)
    return obj


def create_bob(
    vertices: np.ndarray | None = None,
    edges: np.ndarray | None = None,
    faces: np.ndarray | None = None,
    name: str = "NewObject",
    collection: bpy.types.Collection | None = None,
) -> BlenderObject:
    "Create an object but return it wrapped as a BlenderObject"
    return BlenderObject(
        create_object(
            vertices=vertices,
            edges=edges,
            faces=faces,
            name=name,
            collection=collection,
        )
    )
