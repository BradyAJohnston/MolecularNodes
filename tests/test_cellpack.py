import molecularnodes as mn
import pytest
import bpy
from .utils import (
    sample_attribute,
    NumpySnapshotExtension
)
from .constants import (
    data_dir
)

mn.unregister()
mn.register()


@pytest.mark.parametrize('format', ['bcif', 'cif'])
def test_load_cellpack(snapshot_custom: NumpySnapshotExtension, format):
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
    assert snapshot_custom == "\n".join(instance_names)
    assert ens.name == name

    ens.modifiers['MolecularNodes'].node_group.nodes['MN_pack_instances'].inputs['As Points'].default_value = False
    mn.blender.nodes.realize_instances(ens)
    for attribute in ens.data.attributes.keys():
        assert snapshot_custom == sample_attribute(
            ens, attribute, evaluate=True)
    assert snapshot_custom == str(ens['chain_ids'])
