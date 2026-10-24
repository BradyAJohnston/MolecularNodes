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
        mix_shader = s.MixShader(
            factor=mn_color.o.alpha,
            shader=s.TransparentBSDF(),
            shader_001=AmbientOcclusionInternal(
                color=mn_color.o.color, ao_space="Local"
            ),
        )
        with s.Frame("Add Outline"):
            mix_shader_1 = s.MixShader(
                factor=MNFresnel(ior=0.95, factor=1.0),
                shader=mix_shader,
                shader_001=s.Color(),
            )
        _material_output = s.MaterialOutput(mix_shader_1)


MATERIAL = FresnelOutline
MATERIAL_NAME = "Fresnel Outline"

MATERIAL_ASSET_METADATA = {
    "catalog_id": "fc8d3698-34f7-4b7e-8167-a2c0391b171b",
}
