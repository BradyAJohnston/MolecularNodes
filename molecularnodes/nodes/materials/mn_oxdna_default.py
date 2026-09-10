# Material 'MN oxDNA Default', dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# The class is a recipe for the material's shader tree — build recreates the material and runs it into material.node_tree.
from nodebpy import geometry as g
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup


class MNColorInput(CustomShaderGroup):
    _name = "MN Color Input"

    def _build_group(self, tree):
        color = tree.outputs.color("Color", (0.0, 0.0, 0.0, 0.0))
        alpha = tree.outputs.float("Alpha")

        attribute = s.Attribute(attribute_name="Color")
        attribute_1 = s.Attribute(attribute_type="INSTANCER", attribute_name="Color")
        math_1 = g.Math(
            value=attribute_1.o.alpha,
            value_001=0.0,
            operation="GREATER_THAN",
            use_clamp=True,
        )
        attribute.o.alpha * (1.0 - math_1) + attribute_1.o.alpha * math_1 >> alpha
        mix = g.Mix(
            factor_float=math_1,
            a_color=attribute.o.color,
            b_color=attribute_1.o.color,
            data_type="RGBA",
            blend_type="ADD",
            clamp_factor=True,
        )

        mix.o.result_color >> color


class MNOxDNADefault(CustomShaderGroup):
    _name = "Shader Nodetree"

    def _build_group(self, tree):
        group = MNColorInput()
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


MATERIAL = MNOxDNADefault
MATERIAL_NAME = "MN oxDNA Default"
