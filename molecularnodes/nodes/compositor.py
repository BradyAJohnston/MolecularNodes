import bpy
from nodebpy import compositor as c


def default_5x_compositor_node_tree() -> bpy.types.CompositorNodeTree:
    with c.tree("Compositor Nodes") as tree:
        c.RenderLayers().o.image >> tree.outputs.color("Image")

    return tree.tree


def mn_compositor_node_tree() -> bpy.types.CompositorNodeTree:
    with c.tree("MN Compositor") as tree:
        input = tree.inputs.color("Image")
        output = tree.outputs.color("Image")

        image = c.Image(image=bpy.data.images.get("mn_annotations", None))
        c.AlphaOver(input, image, 1.0) >> output

    return tree.tree
