# Material "MN Default.old", dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# The class is a recipe for the material's shader tree — build recreates the material and runs it into material.node_tree.
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup
from ..shader.mn_color import MNColor


class MNDefaultOld(CustomShaderGroup):
    _name = "Shader Nodetree"

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        group = MNColor()
        principled_bsdf = s.PrincipledBSDF(
            base_color=group.o.color,
            alpha=group.o.alpha,
            ior=1.45,
            subsurface_scale=0.05,
            coat_roughness=0.03,
            thin_film_ior=1.33,
            subsurface_method="RANDOM_WALK_LEGACY",
        )
        _material_output = s.MaterialOutput(
            surface=principled_bsdf, is_active_output=True
        )


MATERIAL = MNDefaultOld
MATERIAL_NAME = "MN Default.old"
