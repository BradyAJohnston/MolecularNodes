# Material "Flat Outline", dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# The class is a recipe for the material's shader tree — build recreates the material and runs it into material.node_tree.
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup
from ..shader._shared.mn_fresnel import MNFresnel
from ..shader.mn_color import MNColor


class FlatOutline(CustomShaderGroup):
    _name = "Shader Nodetree"

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        with s.Frame("Shade the Flat Colors with Ambient Occlusion"):
            mix = g.Mix(
                a_color=g.ColorRamp(
                    fac=s.AmbientOcclusion(samples=16).o.ao ** 1.5
                ).o.color,
                b_color=MNColor().o.color,
                data_type="RGBA",
                blend_type="MULTIPLY",
                clamp_factor=True,
            )
        with s.Frame("Stops colors as lights in Cycles"):
            mix_shader = s.MixShader(
                fac=s.LightPath().o.is_camera_ray,
                shader=mix.o.result_color,
                shader_001=s.Emission(color=mix.o.result_color),
            )
        with s.Frame("Add Outline"):
            mix_shader_1 = s.MixShader(
                fac=MNFresnel(ior=0.95, factor=1.0),
                shader=mix_shader,
                shader_001=s.Color(),
            )
        _material_output = s.MaterialOutput(surface=mix_shader_1, is_active_output=True)


MATERIAL = FlatOutline
MATERIAL_NAME = "Flat Outline"

MATERIAL_ASSET_METADATA = {
    "catalog_id": "fc8d3698-34f7-4b7e-8167-a2c0391b171b",
}
