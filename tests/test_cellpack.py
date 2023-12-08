import molecularnodes as mn
import pytest
import bpy
from .utils import (
    sample_attribute_to_string, 
    apply_mods
)
from .constants import (
    test_data_directory
)

@pytest.mark.parametrize('file_format', ['bcif', 'cif'])
def test_cellpack_data(snapshot, file_format):
    object, collection = mn.io.cellpack.parse(
        test_data_directory / f"square1.{file_format}"
    )
    attributes = object.data.attributes.keys()
    for attribute in attributes:
        snapshot.assert_match(
            sample_attribute_to_string(object, attribute), 
            f"att_{attribute}_values.txt"
        )

@pytest.mark.parametrize('file_format', ['bcif', 'cif'])
def test_load_cellpack(snapshot, file_format):
    name = f"Cellpack_{file_format}"
    ens = mn.io.cellpack.load(
        test_data_directory / f"square1.{file_format}", 
        name = name, 
        instance_nodes=False, 
        fraction=0.1
    )
    
    coll = bpy.data.collections[f'cellpack_{name}']
    instance_names = [object.name for object in coll.objects]
    snapshot.assert_match("\n".join(instance_names), "instance_names.txt")
    
    assert ens.name == name
    ens.modifiers['MolecularNodes'].node_group.nodes['MN_pack_instances'].inputs['As Points'].default_value = False
    mn.blender.nodes.realize_instances(ens)
    apply_mods(ens)
    
    for attribute in ens.data.attributes.keys():
        snapshot.assert_match(
            sample_attribute_to_string(ens, attribute), 
            f"att_{attribute}_values.txt"
        )