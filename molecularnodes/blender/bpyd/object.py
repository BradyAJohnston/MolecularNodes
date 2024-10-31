import bpy
import numpy as np
from typing import Optional
from .attribute import (
    evaluate_object,
    AttributeTypes,
    AttributeType,
    Domains,
    DomainType,
)
from . import attribute
from mathutils import Matrix


class ObjectMissingError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


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
    vertices: np.ndarray | None = None,
    edges: np.ndarray | None = None,
    faces: np.ndarray | None = None,
    name: str = "NewObject",
    collection: bpy.types.Collection | None = None,
) -> bpy.types.Object:
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
        bpy.types.Object
            The created object.
    """
    mesh = bpy.data.meshes.new(name)
    if vertices is None:
        vertices = []
    if edges is None:
        edges = []
    if faces is None:
        faces = []
    mesh.from_pydata(vertices=vertices, edges=edges, faces=faces)
    obj = bpy.data.objects.new(name, mesh)
    if not collection:
        collection = bpy.data.collections["Collection"]
    collection.objects.link(obj)
    return obj


def active_object(context: bpy.types.Context = None) -> bpy.types.Object:
    if context is None:
        return bpy.context.active_object

    return context.active_object


class BlenderObject:
    """
    A convenience class for working with Blender objects
    """

    def __init__(self, obj: bpy.types.Object | None):
        if not isinstance(obj, bpy.types.Object):
            raise ValueError(f"{obj} must be a Blender object of type bpy.types.Object")
        self._object = obj

    @property
    def object(self) -> bpy.types.Object:
        obj = self._object
        if obj is None:
            raise ObjectMissingError(
                "Object is deleted and unable to establish a connection with a new Blender Object."
            )
        return obj

    @object.setter
    def object(self, value: bpy.types.Object) -> None:
        self._object = value

    def store_named_attribute(
        self,
        data: np.ndarray,
        name: str,
        atype: str | AttributeType | None = None,
        domain: str | DomainType = Domains.POINT,
    ) -> None:
        """
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
        attribute.store_named_attribute(
            self.object, data=data, name=name, atype=atype, domain=domain
        )
        return self

    def evaluate(self):
        return BlenderObject(evaluate_object(self.object))

    def named_attribute(self, name: str, evaluate: bool = False) -> np.ndarray:
        return attribute.named_attribute(self.object, name=name, evaluate=evaluate)

    def transform_origin(self, matrix: Matrix) -> None:
        self.object.matrix_local = matrix * self.object.matrix_world

    def transform_points(self, matrix: Matrix) -> None:
        self.position = self.position * matrix

    @property
    def selected(self) -> np.ndarray:
        return self.named_attribute(".select_vert")

    @property
    def name(self) -> str:
        obj = self.object
        if obj is None:
            return None

        return obj.name

    @name.setter
    def name(self, value: str) -> None:
        obj = self.object
        if obj is None:
            raise ObjectMissingError
        obj.name = value

    @property
    def position(self) -> np.ndarray:
        return self.named_attribute("position")

    @position.setter
    def position(self, value: np.ndarray) -> None:
        self.store_named_attribute(
            value,
            name="position",
            atype=AttributeTypes.FLOAT_VECTOR,
            domain=Domains.POINT,
        )

    def selected_positions(self, mask: Optional[np.ndarray] = None) -> np.ndarray:
        if mask is not None:
            return self.position[np.logical_and(self.selected, mask)]

        return self.position[self.selected]

    def __len__(self) -> int:
        return len(self.object.data.vertices)
