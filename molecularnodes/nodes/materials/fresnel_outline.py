# Material "Fresnel Outline", dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# The class is a recipe for the material's shader tree — build recreates the material and runs it into material.node_tree.
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup
from ..shader._shared.mn_fresnel import MNFresnel
from ..shader.ambient_occlusion_internal import AmbientOcclusionInternal
from ..shader.mn_color import MNColor


class FresnelOutline(CustomShaderGroup):
    _name = "Shader Nodetree"

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        mn_color = MNColor()
        ambient_occlusion_internal = AmbientOcclusionInternal(color=mn_color.o.color)
        mix_shader = s.MixShader(
            factor=mn_color.o.alpha,
            shader=s.TransparentBSDF(),
            shader_001=ambient_occlusion_internal.o.cycles,
        )
        with s.Frame("Add Outline"):
            color = s.Color()
            mn_fresnel = MNFresnel(ior=0.95, factor=1.0)
            mix_shader_1 = s.MixShader(
                factor=mn_fresnel,
                shader=ambient_occlusion_internal.o.eevee,
                shader_001=color,
            )
            mix_shader_2 = s.MixShader(
                factor=mn_fresnel, shader=mix_shader, shader_001=color
            )
        _material_output = s.MaterialOutput(surface=mix_shader_1, target="EEVEE")
        _material_output_1 = s.MaterialOutput(
            surface=mix_shader_2, is_active_output=True, target="CYCLES"
        )


MATERIAL = FresnelOutline
MATERIAL_NAME = "Fresnel Outline"

MATERIAL_ASSET_METADATA = {
    "catalog_id": "fc8d3698-34f7-4b7e-8167-a2c0391b171b",
}
