import databpy
import numpy as np
import molecularnodes as mn
from molecularnodes.nodes import nodes
from .constants import data_dir


def test_centre_on_selection():
    mol = mn.Molecule.fetch("1cd3", cache=data_dir)

    group = nodes.get_mod(mol.object).node_group = nodes.new_tree()
    link = group.links.new
    node_centre = nodes.add_custom(group, "Centre on Selection")
    node_chain_id = nodes.add_custom(group, "Chain ID", location=[-100, -100])

    link(group.nodes["Group Input"].outputs[0], node_centre.inputs[0])
    link(node_chain_id.outputs["chain_id"], node_centre.inputs["Group ID"])
    link(node_centre.outputs[0], group.nodes["Group Output"].inputs[0])

    old_pos = mol.named_attribute("position", evaluate=False)
    new_pos = mol.named_attribute("position", evaluate=True)

    assert not np.allclose(old_pos, new_pos)

    chain_id = mol.named_attribute("chain_id")
    chain_ids = np.unique(chain_id)

    old_centres = [databpy.centre(old_pos[chain_id == x]) for x in chain_ids]
    new_centres = [databpy.centre(new_pos[chain_id == x]) for x in chain_ids]

    assert not np.allclose(old_centres, new_centres)
    assert np.allclose(
        new_centres, [np.array((0, 0, 0), dtype=float) for x in chain_ids], atol=1e-5
    )


def test_atoms_to_ca_splines(snapshot):
    """
    This test is for if the chains are still connected as continuous, even though they are
    separate by a larger distance
    """
    mol = mn.Molecule.fetch("1HQM")

    group = nodes.get_mod(mol.object).node_group = nodes.new_tree()
    link = group.links.new

    node_atcac = nodes.add_custom(group, "Atoms to CA Curves")
    node_ctm = group.nodes.new("GeometryNodeCurveToMesh")

    link(group.nodes["Group Input"].outputs[0], node_atcac.inputs[0])
    link(node_atcac.outputs[0], node_ctm.inputs[0])
    link(node_ctm.outputs[0], group.nodes["Group Output"].inputs[0])

    pos_pre = mol.named_attribute("position", evaluate=True)

    node_atcac.inputs["Threshold"].default_value = 200

    pos_post = mol.named_attribute("position", evaluate=True)

    assert not np.allclose(pos_pre.shape, pos_post.shape)
