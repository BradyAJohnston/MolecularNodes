import databpy
import numpy as np
import pytest
from molecularnodes.nodes import nodes


def create_debug_group(name="MolecularNodesDebugGroup"):
    group = nodes.new_tree(name=name, fallback=False)
    info = group.nodes.new("GeometryNodeObjectInfo")
    group.links.new(info.outputs["Geometry"], group.nodes["Group Output"].inputs[0])
    return group


custom_selections = [
    ("1, 3, 5-7", np.array((1, 3, 5, 6, 7))),
    ("5, 9-20", np.append(5, np.arange(9, 21))),
    ("1, 7, 8, 9", np.array((1, 7, 8, 9))),
]


@pytest.mark.parametrize("selection", custom_selections)
def test_select_multiple_residues(selection):
    n_atoms = 100
    bob = databpy.create_bob(np.zeros((n_atoms, 3)))
    bob.store_named_attribute(
        data=np.arange(n_atoms) + 1,
        name="res_id",
    )

    mod = nodes.get_mod(bob.object)
    group = nodes.new_tree(fallback=False)
    mod.node_group = group
    sep = group.nodes.new("GeometryNodeSeparateGeometry")
    nodes.insert_last_node(group, sep)

    node_sel_group = nodes.resid_multiple_selection("custom", selection[0])
    node_sel = nodes.add_custom(group, node_sel_group.name)
    group.links.new(node_sel.outputs["Selection"], sep.inputs["Selection"])

    vertices_count = len(bob.evaluate().data.vertices)
    assert vertices_count == len(selection[1])
    assert (bob.named_attribute("res_id", evaluate=True) == selection[1]).all()
