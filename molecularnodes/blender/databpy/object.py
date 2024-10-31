import bpy
import numpy as np
from typing import Union, Optional, Type
from .attribute import (
    store_named_attribute,
    named_attribute,
    Domains,
    DomainType,
)
from mathutils import Matrix


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


class BlenderObject:
    """
    A convenience class for working with Blender objects
    """

    def __init__(self, object: bpy.types.Object):
        self.object = object

    def store_named_attribute(
        self,
        data: np.ndarray,
        name: str,
        dtype: Type | None = None,
        domain: str | DomainType = Domains.POINT,
    ) -> None:
        store_named_attribute(self.object, data=data, name=name, domain=domain)

    def named_attribute(self, name: str, evaluate: bool = False) -> np.ndarray:
        return named_attribute(self.object, name=name, evaluate=evaluate)

    def transform_origin(self, matrix: Matrix) -> None:
        self.object.matrix_local = matrix * self.object.matrix_world

    def transform_points(self, matrix: Matrix) -> None:
        self.position = self.position * matrix

    @property
    def selected(self) -> np.ndarray:
        return named_attribute(self.object, ".select_vert")

    @property
    def position(self) -> np.ndarray:
        return named_attribute(self.object, name="position", evaluate=False)

    @position.setter
    def position(self, value: np.ndarray) -> None:
        store_named_attribute(self.object, "position", value)

    def selected_positions(self, mask: Optional[np.ndarray] = None) -> np.ndarray:
        if mask is not None:
            return self.position[np.logical_and(self.selected, mask)]

        return self.position[self.selected]


def bob(object: Union[bpy.types.Object, BlenderObject]) -> BlenderObject:
    """
    Convenience function to convert a Blender object to a BlenderObject
    """
    if isinstance(object, BlenderObject):
        return object
    elif isinstance(object, bpy.types.Object):
        return BlenderObject(object)
    else:
        raise ValueError(
            f"Unknown object type: {object=}"
            "Expected bpy.types.Object or BlenderObject"
        )
