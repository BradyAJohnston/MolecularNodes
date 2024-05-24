# import bpy
import numpy as np
import pytest
import molecularnodes as mn
from molecularnodes.blender import nodes
import random

from .utils import sample_attribute, NumpySnapshotExtension
from .constants import codes, data_dir
random.seed(6)

mn.unregister()
mn.register()


def test_node_name_format():
    assert mn.blender.nodes.format_node_name(
        "MN_style_cartoon") == "Style Cartoon"
    assert mn.blender.nodes.format_node_name(
        'MN_dna_double_helix'
    ) == 'DNA Double Helix'
    assert mn.blender.nodes.format_node_name(
        'MN_topo_vector_angle'
    ) == 'Topology Vector Angle'


def test_get_nodes():
    bob = mn.io.fetch('4ozs', style='spheres').object

    assert nodes.get_nodes_last_output(bob.modifiers['MolecularNodes'].node_group)[
        0].name == "MN_style_spheres"
    nodes.realize_instances(bob)
    assert nodes.get_nodes_last_output(bob.modifiers['MolecularNodes'].node_group)[
        0].name == "Realize Instances"
    assert nodes.get_style_node(bob).name == "MN_style_spheres"

    bob2 = mn.io.fetch('1cd3', style='cartoon', build_assembly=True).object

    assert nodes.get_nodes_last_output(bob2.modifiers['MolecularNodes'].node_group)[
        0].name == "MN_assembly_1cd3"
    assert nodes.get_style_node(bob2).name == "MN_style_cartoon"


def test_selection():
    chain_ids = [let for let in 'ABCDEFG123456']
    node = nodes.custom_iswitch(
        'test_node', chain_ids, prefix="Chain ", dtype='BOOLEAN')

    input_sockets = nodes.inputs(node)
    for letter, socket in zip(chain_ids, input_sockets.values()):
        assert f"Chain {letter}" == socket.name
        assert socket.default_value is False


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("attribute", ["chain_id", "entity_id"])
def test_selection_working(snapshot_custom: NumpySnapshotExtension, attribute, code):
    mol = mn.io.fetch(code, style='ribbon', cache_dir=data_dir).object
    group = mol.modifiers['MolecularNodes'].node_group
    node_sel = nodes.add_selection(
        group, mol.name, mol[f'{attribute}s'], attribute)

    n = len(node_sel.inputs)

    for i in random.sample(list(range(n)), max(n - 2, 1)):
        node_sel.inputs[i].default_value = True

    nodes.realize_instances(mol)

    assert snapshot_custom == sample_attribute(mol, 'position', evaluate=True)


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("attribute", ["chain_id", 'entity_id'])
def test_color_custom(snapshot_custom: NumpySnapshotExtension, code,  attribute):
    mol = mn.io.fetch(code, style='ribbon', cache_dir=data_dir).object

    group_col = mn.blender.nodes.custom_iswitch(
        name=f'MN_color_entity_{mol.name}',
        iter_list=mol[f'{attribute}s'],
        field=attribute,
        dtype='RGBA'
    )
    group = mol.modifiers['MolecularNodes'].node_group
    node_col = mn.blender.nodes.add_custom(
        group, group_col.name, [0, -200])
    group.links.new(
        node_col.outputs[0],
        group.nodes['MN_color_set'].inputs['Color']
    )

    assert snapshot_custom == sample_attribute(mol, 'Color', n=50)


def test_custom_resid_selection():
    node = mn.blender.nodes.resid_multiple_selection(
        'new_node', '1, 5, 10-20, 40-100')
    numbers = [1, 5, 10, 20, 40, 100]
    assert len(nodes.outputs(node)) == 2
    counter = 0
    for item in node.interface.items_tree:
        if item.in_out == "INPUT":
            assert item.default_value == numbers[counter]
            counter += 1


def test_op_custom_color():
    mol = mn.io.load(data_dir / '1cd3.cif').object
    mol.select_set(True)
    group = mn.blender.nodes.custom_iswitch(
        name=f'MN_color_chain_{mol.name}',
        iter_list=mol['chain_ids'],
        dtype='RGBA'
    )

    assert group
    assert group.interface.items_tree['G'].name == 'G'
    assert group.interface.items_tree[-1].name == 'G'
    assert group.interface.items_tree[0].name == 'Color'


def test_color_lookup_supplied():
    col = mn.color.random_rgb(6)
    name = 'test'
    node = mn.blender.nodes.custom_iswitch(
        name=name,
        iter_list=range(10, 20),
        dtype='RGBA',
        default_values=[col for i in range(10)],
        start=10
    )
    assert node.name == name
    for item in nodes.inputs(node).values():
        assert np.allclose(np.array(item.default_value), col)

    node = mn.blender.nodes.custom_iswitch(
        name='test2',
        iter_list=range(10, 20),
        dtype='RGBA',
        start=10
    )
    for item in nodes.inputs(node).values():
        assert not np.allclose(np.array(item.default_value), col)


def test_color_chain(snapshot_custom: NumpySnapshotExtension):
    mol = mn.io.load(data_dir / '1cd3.cif', style='cartoon').object
    group_col = mn.blender.nodes.custom_iswitch(
        name=f'MN_color_chain_{mol.name}',
        iter_list=mol['chain_ids'],
        dtype='RGBA'
    )
    group = mol.modifiers['MolecularNodes'].node_group
    node_col = mn.blender.nodes.add_custom(group, group_col.name, [0, -200])
    group.links.new(node_col.outputs[0],
                    group.nodes['MN_color_set'].inputs['Color'])

    assert snapshot_custom == sample_attribute(mol, 'Color')


def test_color_entity(snapshot_custom: NumpySnapshotExtension):
    mol = mn.io.fetch('1cd3', style='cartoon').object
    group_col = mn.blender.nodes.custom_iswitch(
        name=f'MN_color_entity_{mol.name}',
        iter_list=mol['entity_ids'],
        dtype='RGBA',
        field='entity_id'
    )
    group = mol.modifiers['MolecularNodes'].node_group
    node_col = mn.blender.nodes.add_custom(group, group_col.name, [0, -200])
    group.links.new(node_col.outputs[0],
                    group.nodes['MN_color_set'].inputs['Color'])

    assert snapshot_custom == sample_attribute(mol, 'Color')


def get_links(sockets):
    links = []
    for socket in sockets:
        for link in socket.links:
            yield link


def test_change_style():
    model = mn.io.fetch('1cd3', style='cartoon').object
    style_node_1 = nodes.get_style_node(model).name
    mn.blender.nodes.change_style_node(model, 'ribbon')
    style_node_2 = nodes.get_style_node(model).name

    assert style_node_1 != style_node_2

    for style in ['ribbon', 'cartoon', 'presets', 'ball_and_stick', 'surface']:
        style_node_1 = nodes.get_style_node(model)
        links_in_1 = [link.from_socket.name for link in get_links(
            style_node_1.inputs)]
        links_out_1 = [link.from_socket.name for link in get_links(
            style_node_1.outputs)]

        nodes.change_style_node(model, style)
        style_node_2 = nodes.get_style_node(model)
        links_in_2 = [link.from_socket.name for link in get_links(
            style_node_2.inputs)]
        links_out_2 = [link.from_socket.name for link in get_links(
            style_node_2.outputs)]

        assert len(links_in_1) == len(links_in_2)
        assert len(links_out_1) == len(links_out_2)


def test_node_topology(snapshot_custom: NumpySnapshotExtension):
    mol = mn.io.fetch('1bna', del_solvent=False).object

    group = nodes.get_mod(mol).node_group

    group.links.new(group.nodes['Group Input'].outputs[0],
                    group.nodes['Group Output'].inputs[0])
    node_att = group.nodes.new('GeometryNodeStoreNamedAttribute')
    node_att.inputs[2].default_value = 'test_attribute'
    nodes.insert_last_node(group, node_att)
    node_names = [
        node['name']
        for node in mn.ui.node_info.menu_items['topology']
        if not node == "break"
    ]
    for node_name in node_names:
        # exclude these particular nodes, as they aren't field nodes and so we shouldn't
        # be testing them here. Will create their own particular tests later
        if 'backbone' in node_name or 'bonds' in node_name:
            continue
        node_topo = nodes.add_custom(group, node_name,
                                     location=[x - 300 for x in node_att.location])

        if node_name == "MN_topo_point_mask":
            node_topo.inputs['atom_name'].default_value = 61

        type_to_data_type = {
            'VECTOR': 'FLOAT_VECTOR',
            'VALUE': 'FLOAT',
            'BOOLEAN': 'BOOLEAN',
            'INT': 'INT',
            'RGBA': 'FLOAT_COLOR',
            'ROTATION': 'QUATERNION'
        }

        for output in node_topo.outputs:
            node_att.data_type = type_to_data_type[output.type]
            input = node_att.inputs['Value']

            for link in input.links:
                group.links.remove(link)

            group.links.new(output, input)

            assert snapshot_custom == mn.blender.obj.get_attribute(
                mol, 'test_attribute', evaluate=True
            )


def test_compute_backbone(snapshot_custom: NumpySnapshotExtension):
    mol = mn.io.fetch('1CCN', del_solvent=False).object

    group = nodes.get_mod(mol).node_group

    group.links.new(group.nodes['Group Input'].outputs[0],
                    group.nodes['Group Output'].inputs[0])
    node_att = group.nodes.new('GeometryNodeStoreNamedAttribute')
    node_att.inputs[2].default_value = 'test_attribute'
    node_backbone = nodes.add_custom(group, 'MN_topo_compute_backbone')
    nodes.insert_last_node(group, node_backbone)
    nodes.insert_last_node(group, node_att)
    node_names = ['MN_topo_backbone']
    for node_name in node_names:
        node_topo = nodes.add_custom(group, node_name,
                                     location=[x - 300 for x in node_att.location])

        if node_name == "MN_topo_point_mask":
            node_topo.inputs['atom_name'].default_value = 61

        type_to_data_type = {
            'VECTOR': 'FLOAT_VECTOR',
            'VALUE': 'FLOAT',
            'BOOLEAN': 'BOOLEAN',
            'INT': 'INT',
            'RGBA': 'FLOAT_COLOR',
            'ROTATION': 'QUATERNION'
        }

        for output in node_topo.outputs:
            node_att.data_type = type_to_data_type[output.type]
            input = node_att.inputs['Value']

            for link in input.links:
                group.links.remove(link)

            group.links.new(output, input)

            assert snapshot_custom == mn.blender.obj.get_attribute(
                mol, 'test_attribute', evaluate=True
            )

        for angle in ['Phi', 'Psi']:
            output = node_backbone.outputs[angle]
            node_att.data_type = type_to_data_type[output.type]
            input = node_att.inputs['Value']

            for link in input.links:
                group.links.remove(link)

            group.links.new(output, input)

            assert snapshot_custom == mn.blender.obj.get_attribute(
                mol, 'test_attribute', evaluate=True
            )


def test_topo_bonds():
    mol = mn.io.fetch('1BNA', del_solvent=True, style=None).object
    group = nodes.get_mod(mol).node_group = nodes.new_group()

    # add the node that will break bonds, set the cutoff to 0
    node_break = nodes.add_custom(group, 'MN_topo_bonds_break')
    nodes.insert_last_node(group, node=node_break)
    node_break.inputs['Cutoff'].default_value = 0

    # compare the number of edges before and after deleting them with
    bonds = mol.data.edges
    no_bonds = mn.blender.obj.evaluated(mol).data.edges
    assert len(bonds) > len(no_bonds)
    assert len(no_bonds) == 0

    # add the node to find the bonds, and ensure the number of bonds pre and post the nodes
    # are the same (other attributes will be different, but for now this is good)
    node_find = nodes.add_custom(group, 'MN_topo_bonds_find')
    nodes.insert_last_node(group, node=node_find)
    bonds_new = mn.blender.obj.evaluated(mol).data.edges
    assert len(bonds) == len(bonds_new)
