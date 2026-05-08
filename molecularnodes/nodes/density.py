import bpy
from nodebpy import geometry as g
from molecularnodes import material


def style_density_iso():
    with g.tree("Style Density ISO Surface") as tree:
        volume = tree.inputs.geometry("Volume")
        visible = tree.inputs.boolean("Visible")
        smooth = tree.inputs.boolean("Smooth")
        iso = tree.inputs.float("ISO Value")
        col_pos = tree.inputs.color("Positive Color")
        col_neg = tree.inputs.color("Negative Color")
        mat = tree.inputs.material("Material")

        volume = visible.switch.geometry(None, volume)
        geom = g.JoinGeometry(
            [
                g.VolumeToMesh(volume, voxel_size=val)
                >> g.StoreNamedAttribute.point.color(name="Color", value=col)
                >> g.SetMaterial(material=mat)
                for val, col in ((iso, col_pos), (-iso, col_neg))
            ]
        ) >>
        g.SetShadeSmooth.face(shade_smooth=smooth)

        # del_sel =
        #


        geom >> tree.outputs.geometry("Geometry")
