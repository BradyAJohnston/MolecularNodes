import bpy
from nodebpy import compositor as c


def compositor_node_tree() -> bpy.types.CompositorNodeTree:
    with c.tree("MN Compositor") as tree:
        (
            tree.inputs.color("Image")
            >> c.AlphaOver(
                background=...,
                foreground=c.Image(bpy.data.images.get("mn_annotations", None)),
            )
            >> tree.outputs.color("Image")
        )

    return tree.tree
