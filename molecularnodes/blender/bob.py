from typing import Optional, Union
import biotite.structure as struc
import bpy
import numpy as np
from mathutils import Matrix

from . import mesh

# computes the RMSD and moves the current object
# allows for usage as such:
# o = BlenderObject(bpy.data.objects['mobile_object'])
# o.align_to(bpy.data.objects['static_object'])


class BlenderObject:
    """
    A convenience class for working with Blender objects
    """

    def __init__(self, object: bpy.types.Object):
        self.object = object

    def store_named_attribute(
        self, name: str, data: np.ndarray, domain: str = "POINT"
    ) -> None:
        mesh.store_named_attribute(self.object, name=name, data=data, domain=domain)

    def named_attribute(self, name, evaluate=False) -> np.ndarray:
        return mesh.named_attribute(self.object, name=name, evaluate=evaluate)

    def transform_origin(self, matrix: Matrix) -> None:
        self.object.matrix_local = matrix * self.object.matrix_world
    
    def transform_points(self, matrix: Matrix) -> None:
        self.position = self.

    @property
    def selected(self) -> np.ndarray:
        return mesh.named_attribute(self.object, ".select_vert")

    @property
    def position(self) -> np.ndarray:
        return mesh.named_attribute(self.object, name="position", evaluate=False)

    @position.setter
    def position(self, value: np.ndarray) -> None:
        mesh.store_named_attribute(self.object, "position", value)

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
            f"Unknown object type: {object}"
            "Expected bpy.types.Object or BlenderObject"
        )


def align_rmsd(
    fixed: Union[bpy.types.Object, BlenderObject],
    mobile: Union[bpy.types.Object, BlenderObject],
    mask_fixed: Optional[np.ndarray] = None,
    mask_mobile: Optional[np.ndarray] = None,
    transform = "origin"
) -> None:
    fixed = bob(fixed)
    mobile = bob(mobile)

    pos_fixed = fixed.selected_positions(mask_fixed)
    pos_mobile = mobile.selected_positions(mask_mobile)

    if len(pos_fixed) != len(pos_mobile):
        if len(pos_mobile) > len(pos_fixed):
            pos_mobile = pos_mobile[0 : len(pos_fixed)]
        else:
            pos_fixed = pos_fixed[0 : len(pos_mobile)]

    fitted, transformation = struc.superimpose(fixed=pos_fixed, mobile=pos_mobile)

    if transform == "origin":
        mobile.transform_origin(Matrix(transformation.as_matrix()))
    elif transform == "points":
        mobile.transform_points(Matrix(transformation.as_matrix()))

    # apply the transformation to the current object's positions
    self.position = transformation.apply(self.position)
