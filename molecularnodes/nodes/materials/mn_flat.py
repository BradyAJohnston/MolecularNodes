# Material "MN Flat", dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# The class is a recipe for the material's shader tree — build recreates the material and runs it into material.node_tree.
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup
from ..shader.flat import Flat


class MNFlat(CustomShaderGroup):
    _name = "Shader Nodetree"

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        _material_output = Flat() >> s.MaterialOutput(is_active_output=True)


MATERIAL = MNFlat
MATERIAL_NAME = "MN Flat"

MATERIAL_ASSET_METADATA = {
    "catalog_id": "fc8d3698-34f7-4b7e-8167-a2c0391b171b",
}
