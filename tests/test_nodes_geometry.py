import numpy as np

import molecularnodes as mn
from molecularnodes.blender import mesh, nodes

from .constants import data_dir


def test_centre_on_selection():
    mol = mn.entities.fetch("1cd3", cache_dir=data_dir, style=None)

    group = nodes.get_mod(mol.object).node_group = nodes.new_group()
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

    old_centres = [mesh.centre(old_pos[chain_id == x]) for x in chain_ids]
    new_centres = [mesh.centre(new_pos[chain_id == x]) for x in chain_ids]

    assert not np.allclose(old_centres, new_centres)
    assert np.allclose(
        new_centres, [np.array((0, 0, 0), dtype=float) for x in chain_ids], atol=1e-5
    )
