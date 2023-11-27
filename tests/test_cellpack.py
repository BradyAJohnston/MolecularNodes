import molecularnodes as mn
import bpy
from .utils import (
    sample_attribute_to_string, 
    apply_mods
)
from .constants import (
    test_data_directory
)

def test_cellpack_data(snapshot):
    object, collection = mn.io.cellpack.parse(
        test_data_directory / "synvesicle_2-no_bonds.bcif"
    )
    attributes = object.data.attributes.keys()
    for attribute in attributes:
        snapshot.assert_match(
            sample_attribute_to_string(object, attribute), 
            f"att_{attribute}_values.txt"
        )

def test_load_cellpack(snapshot):
    name = "Cellpack"
    mn.io.cellpack.load(
        test_data_directory / "synvesicle_2-no_bonds.bcif", 
        name = name, 
        instance_nodes=False, 
        fraction=0.1
    )
    
    obj = bpy.data.objects[name]
    coll = bpy.data.collections[f'cellpack_{name}']
    instance_names = [object.name for object in coll.objects]
    snapshot.assert_match("\n".join(instance_names), "instance_names.txt")
    
    assert obj.name == name
    obj.modifiers['MolecularNodes'].node_group.nodes['MN_pack_instances'].inputs['As Points'].default_value = False
    mn.blender.nodes.realize_instances(obj)
    apply_mods(obj)
    
    for attribute in obj.data.attributes.keys():
        snapshot.assert_match(
            sample_attribute_to_string(obj, attribute), 
            f"att_{attribute}_values.txt"
        )