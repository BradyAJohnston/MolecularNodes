import databpy
import numpy as np
import molecularnodes as mn
from molecularnodes.nodes import nodes
from .constants import data_dir


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
