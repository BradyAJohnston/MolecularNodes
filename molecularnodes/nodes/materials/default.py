# Material "Default", dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# The class is a recipe for the material's shader tree — build recreates the material and runs it into material.node_tree.
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup
from ..shader.color_ao import ColorAO
from ..shader.mn_color import MNColor


class Default(CustomShaderGroup):
    _name = "Shader Nodetree"

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        group = MNColor()
        principled_bsdf = s.PrincipledBSDF(
            base_color=ColorAO(color=group.o.color, distance=0.5, exponent=1.0),
            alpha=group.o.alpha,
            roughness=0.2636364,
            ior=1.45,
            subsurface_scale=0.0,
            coat_roughness=0.03,
            thin_film_ior=1.33,
            subsurface_method="RANDOM_WALK_LEGACY",
        )
        _material_output = s.MaterialOutput(
            surface=principled_bsdf, is_active_output=True
        )


MATERIAL = Default
MATERIAL_NAME = "Default"

MATERIAL_ASSET_METADATA = {
    "catalog_id": "fc8d3698-34f7-4b7e-8167-a2c0391b171b",
}
