# Material "MN Ambient Occlusion", dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# The class is a recipe for the material's shader tree — build recreates the material and runs it into material.node_tree.
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup
from ..shader.color_ao import ColorAO
from ..shader.mn_color import MNColor


class MNAmbientOcclusion(CustomShaderGroup):
    _name = "Shader Nodetree"

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        _material_output = s.MaterialOutput(
            surface=s.Emission(color=ColorAO(color=MNColor().o.color)),
            is_active_output=True,
        )


MATERIAL = MNAmbientOcclusion
MATERIAL_NAME = "MN Ambient Occlusion"

MATERIAL_PROPERTIES = {
    "surface_render_method": "BLENDED",
    "use_transparency_overlap": False,
    "blend_method": "BLEND",
}
