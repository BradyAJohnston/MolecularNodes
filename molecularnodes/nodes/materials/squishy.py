# Material "Squishy", dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# The class is a recipe for the material's shader tree — build recreates the material and runs it into material.node_tree.
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup
from ..shader.mn_color import MNColor


class Squishy(CustomShaderGroup):
    _name = "Shader Nodetree"

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        group = MNColor()
        principled_bsdf = s.PrincipledBSDF(
            base_color=group.o.color,
            alpha=group.o.alpha,
            roughness=1.0,
            ior=1.05,
            diffuse_roughness=1.0,
            subsurface_weight=1.0,
            subsurface_scale=0.2,
            coat_weight=1.0,
            coat_roughness=0.24545455,
            thin_film_ior=1.33,
        )
        _material_output = s.MaterialOutput(
            surface=principled_bsdf, is_active_output=True
        )


MATERIAL = Squishy
MATERIAL_NAME = "Squishy"

MATERIAL_ASSET_METADATA = {
    "catalog_id": "fc8d3698-34f7-4b7e-8167-a2c0391b171b",
}
