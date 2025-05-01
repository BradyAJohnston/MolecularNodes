import bpy
import numpy as np
import molecularnodes as mn
from molecularnodes.nodes.nodes import get_style_node
from molecularnodes.nodes.styles import (
    StyleBallandStick,
    StyleCartoon,
    StyleRibbon,
    StyleSpheres,
    StyleSticks,
    StyleSurface,
)


DONT_COMPARE = {"Atoms", "Selection", "Material"}


def assess_node_equivalency(name, style):
    # get the names and default values of the blender style ndoes
    mol = mn.Molecule.fetch("4ozs").add_style(name)

    style_node = get_style_node(mol.object)
    print(name)
    for input in style_node.inputs:
        print(input.name, input.type)

    blender_inputs = [
        [input.name, input.default_value]
        for input in style_node.inputs
            if input.type != "GEOMETRY" and
            input.name not in DONT_COMPARE

    ]
    blender_names = set(name for [name, _] in blender_inputs)

    # get the style class name
    style_class = style()
    style_class_bnames = set(sc.get("blendername") for sc in style_class.portdata)

    # check names bidirectionally
    for bname in blender_names:
        assert bname in style_class_bnames, (
            f"MN Style {name} has field {bname} that is not found in the styles class"
        )
    for sname in style_class_bnames:
        assert sname in blender_names, (
            f"Internal Style class for {name} has field {bname} that is not found in the upstream MN Style Node"
        )

    for [bname, bvalue] in blender_inputs:
        for pdata in style_class.portdata:
            if pdata.get("blendername") == bname:
                local_name = pdata.get("name")
                local_val = getattr(style_class, local_name)

                # floats come from C++ and are artificaially long
                if isinstance(bvalue, float):
                    bvalue = round(bvalue, 2)

                assert local_val == bvalue, (
                    "Values {bvalue} and {local_val} are not equivalent"
                )

    # equivalency check.
    assert style_class_bnames.difference(blender_names) == set()
    assert blender_names.difference(style_class_bnames) == set()

    bpy.data.objects.remove(mol.object)


def test_styles():
    assess_node_equivalency("ball+stick", StyleBallandStick)
    assess_node_equivalency("cartoon", StyleCartoon)
    assess_node_equivalency("ribbon", StyleRibbon)
    assess_node_equivalency("spheres", StyleSpheres)
    assess_node_equivalency("sticks", StyleSticks)
    assess_node_equivalency("surface", StyleSurface)
