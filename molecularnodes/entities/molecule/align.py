import biotite.structure as struc
import bpy
import numpy as np

from ...blender import mesh

# computes the RMSD and moves the current object
# allows for usage as such:
# o = BlenderObject(bpy.data.objects['mobile_object'])
# o.align_to(bpy.data.objects['static_object'])


class BlenderObject:
    def __init__(self, object: bpy.types.Object):
        self.object = object

    def store_named_attribute(
        self, name: str, data: np.ndarray, domain: str = "POINT"
    ) -> None:
        mesh.store_named_attribute(self.object, name=name, data=data, domain=domain)

    def named_attribute(self, name, evaluate=False) -> np.ndarray:
        return mesh.named_attribute(self.object, name=name, evaluate=evaluate)

    @property
    def selected(self) -> np.ndarray:
        return mesh.named_attribute(self.object, ".select_vert")

    @property
    def position(self) -> np.ndarray:
        return mesh.named_attribute(self.object)

    @position.setter
    def position(self, value: np.ndarray) -> None:
        mesh.store_named_attribute(self.object, "position", value)

    def selected_positions(self, ca_only=False) -> np.ndarray:
        return self.position[
            np.logical_and(self.selected, self.named_attribute("is_alpha_carbon"))
        ]

    def align_to(self, fixed: bpy.types.Object, ca_only=True) -> None:
        fixed = BlenderObject(fixed)

        pos_fixed = fixed.selected_positions(ca_only=ca_only)
        pos_mobile = self.selected_positions(ca_only=ca_only)

        if len(pos_fixed) != len(pos_mobile):
            if len(pos_mobile) > len(pos_fixed):
                pos_mobile = pos_mobile[0 : len(pos_fixed)]
            else:
                pos_fixed = pos_fixed[0 : len(pos_mobile)]

        fitted, transformation = struc.superimpose(fixed=pos_fixed, mobile=pos_mobile)

        # apply the transformation to the current object's positions
        self.position = transformation.apply(self.position)
