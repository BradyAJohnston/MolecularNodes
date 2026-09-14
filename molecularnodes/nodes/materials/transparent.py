# Material "Transparent", dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# The class is a recipe for the material's shader tree — build recreates the material and runs it into material.node_tree.
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup
from ..shader.transparent_outline_internal import TransparentOutlineInternal


class Transparent(CustomShaderGroup):
    _name = "Shader Nodetree"

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        _material_output = TransparentOutlineInternal(
            transparency=0.8083333,
            menu="Fresnel",
            outline_color=(0.0, 0.0, 0.0, 1.0),
            threshold=6.2,
        ) >> s.MaterialOutput(is_active_output=True)


MATERIAL = Transparent
MATERIAL_NAME = "Transparent"

MATERIAL_PROPERTIES = {
    "surface_render_method": "BLENDED",
    "use_backface_culling": True,
    "use_transparency_overlap": False,
    "blend_method": "BLEND",
}

MATERIAL_ASSET_METADATA = {
    "catalog_id": "fc8d3698-34f7-4b7e-8167-a2c0391b171b",
}
