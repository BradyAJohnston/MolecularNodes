import bpy
import numpy as np


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


def evaluate_object(obj):
    "Return an object which has the modifiers evaluated."
    obj.update_tag()
    return obj.evaluated_get(bpy.context.evaluated_depsgraph_get())

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
