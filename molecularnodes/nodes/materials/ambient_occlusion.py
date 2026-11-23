# Material "Ambient Occlusion", dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# The class is a recipe for the material's shader tree — build recreates the material and runs it into material.node_tree.
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup
from ..shader.ambient_occlusion_internal import AmbientOcclusionInternal
from ..shader.mn_color import MNColor


class AmbientOcclusion(CustomShaderGroup):
    _name = "Shader Nodetree"

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        mn_color = MNColor()
        ambient_occlusion_internal = AmbientOcclusionInternal(color=mn_color.o.color)
        _material_output = s.MaterialOutput.eevee(ambient_occlusion_internal.o.eevee)
        mix_shader = s.MixShader(
            factor=mn_color.o.alpha,
            shader=s.TransparentBSDF(),
            shader_001=ambient_occlusion_internal.o.cycles,
        )
        _material_output_1 = s.MaterialOutput.cycles(mix_shader, is_active_output=True)


MATERIAL = AmbientOcclusion
MATERIAL_NAME = "Ambient Occlusion"

MATERIAL_PROPERTIES = {
    "surface_render_method": "BLENDED",
    "use_transparency_overlap": False,
    "blend_method": "BLEND",
}

MATERIAL_ASSET_METADATA = {
    "catalog_id": "fc8d3698-34f7-4b7e-8167-a2c0391b171b",
}
