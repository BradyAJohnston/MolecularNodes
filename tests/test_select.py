import molecularnodes as mn
from molecularnodes.blender import nodes
import bpy
import numpy as np
import pytest


def create_debug_group(name='MolecularNodesDebugGroup'):
    group = nodes.new_group(name=name, fallback=False)
    info = group.nodes.new('GeometryNodeObjectInfo')
    group.links.new(info.outputs['Geometry'],
                    group.nodes['Group Output'].inputs[0])
    return group


def evaluate(object):
    object.update_tag()
    dg = bpy.context.evaluated_depsgraph_get()
    return object.evaluated_get(dg)


custom_selections = [
    ('1, 3, 5-7', np.array((1, 3, 5, 6, 7))),
    ('5, 9-20', np.append(5, np.arange(9, 21))),
    ('1, 7, 8, 9', np.array((1, 7, 8, 9)))
]


@pytest.mark.parametrize('selection', custom_selections)
def test_select_multiple_residues(selection):
    n_atoms = 100
    object = mn.blender.obj.create_object(np.zeros((n_atoms, 3)))
    mn.blender.obj.set_attribute(object, 'res_id', np.arange(n_atoms) + 1)

    mod = nodes.get_mod(object)
    group = nodes.new_group(fallback=False)
    mod.node_group = group
    sep = group.nodes.new('GeometryNodeSeparateGeometry')
    nodes.insert_last_node(group, sep)

    node_sel_group = nodes.resid_multiple_selection('custom', selection[0])
    node_sel = nodes.add_custom(group, node_sel_group.name)
    group.links.new(node_sel.outputs['Selection'], sep.inputs['Selection'])

    vertices_count = len(mn.blender.obj.evaluated(object).data.vertices)
    assert vertices_count == len(selection[1])
    assert (mn.blender.obj.get_attribute(
        mn.blender.obj.evaluated(object), 'res_id') == selection[1]).all()
