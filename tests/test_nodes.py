import bpy
import numpy as np
import pytest
import molecularnodes as mn
from molecularnodes.blender import nodes
import random
import tempfile

from . import utils
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
    mol = mn.io.fetch('4ozs', style='spheres').object

    assert nodes.get_nodes_last_output(mol.modifiers['MolecularNodes'].node_group)[
        0].name == "MN_style_spheres"
    nodes.realize_instances(mol)
    assert nodes.get_nodes_last_output(mol.modifiers['MolecularNodes'].node_group)[
        0].name == "Realize Instances"
    assert nodes.get_style_node(mol).name == "MN_style_spheres"

    mol2 = mn.io.fetch('1cd3', style='cartoon', build_assembly=True).object

    assert nodes.get_nodes_last_output(mol2.modifiers['MolecularNodes'].node_group)[
        0].name == "MN_assembly_1cd3"
    assert nodes.get_style_node(mol2).name == "MN_style_cartoon"


def test_selection():
    chain_ids = [let for let in 'ABCDEFG123456']
    node = nodes.chain_selection('test_node', chain_ids, label_prefix="Chain ")

    input_sockets = nodes.inputs(node)
    for letter, socket in zip(chain_ids, input_sockets.values()):
        assert f"Chain {letter}" == socket.name
        assert socket.default_value is False


with tempfile.TemporaryDirectory() as temp:
    @pytest.mark.parametrize("code", codes)
    @pytest.mark.parametrize("format", ['mmtf', 'cif'])
    @pytest.mark.parametrize("attribute", ["chain_id", "entity_id"])
    def test_selection_working(snapshot, attribute, code, format):
        mol = mn.io.fetch(code, style='ribbon',
                          cache_dir=temp, format=format).object
        group = mol.modifiers['MolecularNodes'].node_group
        node_sel = nodes.add_selection(
            group, mol.name, mol[f'{attribute}s'], attribute)

        n = len(node_sel.inputs)

        for i in random.sample(list(range(n)), max(n - 2, 1)):
            node_sel.inputs[i].default_value = True

        nodes.realize_instances(mol)

        snapshot.assert_match(
            utils.sample_attribute_to_string(utils.evaluate(mol), 'position'),
            "position.txt"
        )

    @pytest.mark.parametrize("code", codes)
    @pytest.mark.parametrize("format", ['mmtf', 'cif'])
    @pytest.mark.parametrize("attribute", ["chain_id", "entity_id"])
    def test_color_custom(snapshot, code, format, attribute):
        mol = mn.io.fetch(code, style='ribbon',
                          format=format, cache_dir=temp).object

        group_col = mn.blender.nodes.chain_color(
            f'MN_color_entity_{mol.name}', input_list=mol[f'{attribute}s'], field=attribute)
        group = mol.modifiers['MolecularNodes'].node_group
        node_col = mn.blender.nodes.add_custom(
            group, group_col.name, [0, -200])
        group.links.new(
            node_col.outputs[0], group.nodes['MN_color_set'].inputs['Color'])

        snapshot.assert_match(
            utils.sample_attribute_to_string(
                utils.evaluate(mol), 'Color', n=500),
            'color.txt'
        )


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
    group = mn.blender.nodes.chain_color(
        f'MN_color_chain_{mol.name}', input_list=mol['chain_ids'])

    assert group
    assert group.interface.items_tree['Chain G'].name == 'Chain G'
    assert group.interface.items_tree[-1].name == 'Chain G'
    assert group.interface.items_tree[0].name == 'Color'


def test_color_chain(snapshot):
    mol = mn.io.load(data_dir / '1cd3.cif', style='cartoon').object
    group_col = mn.blender.nodes.chain_color(
        f'MN_color_chain_{mol.name}', input_list=mol['chain_ids'])
    group = mol.modifiers['MolecularNodes'].node_group
    node_col = mn.blender.nodes.add_custom(group, group_col.name, [0, -200])
    group.links.new(node_col.outputs[0],
                    group.nodes['MN_color_set'].inputs['Color'])

    snapshot.assert_match(
        utils.sample_attribute_to_string(mol, 'Color', n=500),
        'color_chain_values.txt'
    )


def test_color_entity(snapshot):
    mol = mn.io.fetch('1cd3', style='cartoon').object
    group_col = mn.blender.nodes.chain_color(
        f'MN_color_entity_{mol.name}', input_list=mol['entity_ids'], field='entity_id')
    group = mol.modifiers['MolecularNodes'].node_group
    node_col = mn.blender.nodes.add_custom(group, group_col.name, [0, -200])
    group.links.new(node_col.outputs[0],
                    group.nodes['MN_color_set'].inputs['Color'])

    snapshot.assert_match(
        utils.sample_attribute_to_string(mol, 'Color', n=500),
        'color_entity_values.txt'
    )


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


def test_node_topology(snapshot):
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
        if 'backbone' in node_name:
            continue
        node_topo = nodes.add_custom(group, node_name,
                                     location=[x - 300 for x in node_att.location])

        if node_name == "MN_topo_point_mask":
            node_topo.inputs['atom_name'].default_value = 61

        dtype_pos_lookup = {
            'VECTOR': 3,  # Value_Vector
            'VALUE': 4,  # Value_Float
            'RGBA': 5,  # 'Value_Color',
            'BOOLEAN': 6,  # 'Value_Bool',
            'INT': 7,  # 'Value_Int',
            'ROTATION': 8,  # 'Value_Rotation'
        }

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
            input = node_att.inputs[dtype_pos_lookup[output.type]]

            for link in input.links:
                group.links.remove(link)

            group.links.new(output, input)

            snapshot.assert_match(
                np.array2string(
                    mn.blender.obj.get_attribute(
                        utils.evaluate(mol), 'test_attribute'),
                    threshold=10000,
                    precision=3
                ),
                f'{node_name}_{output.name}.txt'

            )


def test_compute_backbone(snapshot):
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

        dtype_pos_lookup = {
            'VECTOR': 3,  # Value_Vector
            'VALUE': 4,  # Value_Float
            'RGBA': 5,  # 'Value_Color',
            'BOOLEAN': 6,  # 'Value_Bool',
            'INT': 7,  # 'Value_Int',
            'ROTATION': 8,  # 'Value_Rotation'
        }

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
            input = node_att.inputs[dtype_pos_lookup[output.type]]

            for link in input.links:
                group.links.remove(link)

            group.links.new(output, input)

            snapshot.assert_match(
                np.array2string(
                    mn.blender.obj.get_attribute(
                        utils.evaluate(mol), 'test_attribute'),
                    threshold=10000,
                    precision=3
                ),
                f'{node_name}_{output.name}.txt'

            )

        for angle in ['Phi', 'Psi']:
            output = node_backbone.outputs[angle]
            node_att.data_type = type_to_data_type[output.type]
            input = node_att.inputs[dtype_pos_lookup[output.type]]

            for link in input.links:
                group.links.remove(link)

            group.links.new(output, input)

            snapshot.assert_match(
                np.array2string(
                    mn.blender.obj.get_attribute(
                        utils.evaluate(mol), 'test_attribute'),
                    threshold=10000,
                    precision=3
                ),
                f'{node_name}_{output.name}.txt'

            )
