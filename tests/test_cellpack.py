import molecularnodes as mn
import pytest
import bpy
from .utils import (
    sample_attribute_to_string
)
from .constants import (
    data_dir
)

mn.unregister()
mn.register()


@pytest.mark.parametrize('format', ['bcif', 'cif'])
def test_load_cellpack(snapshot, format):
    bpy.ops.wm.read_homefile(app_template="")
    name = f"Cellpack_{format}"
    ens = mn.io.cellpack.load(
        data_dir / f"square1.{format}",
        name=name,
        node_setup=False,
        fraction=0.1
    )

    coll = bpy.data.collections[f'cellpack_{name}']
    instance_names = [object.name for object in coll.objects]
    snapshot.assert_match("\n".join(instance_names), "instance_names.txt")
    assert ens.name == name

    ens.modifiers['MolecularNodes'].node_group.nodes['MN_pack_instances'].inputs['As Points'].default_value = False
    mn.blender.nodes.realize_instances(ens)
    for attribute in ens.data.attributes.keys():
        snapshot.assert_match(
            sample_attribute_to_string(ens, attribute, evaluate=True),
            f"att_{attribute}_values.txt"
        )
    snapshot.assert_match(str(ens['chain_ids']), "chain_ids.txt")
