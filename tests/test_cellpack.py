import molecularnodes as mn
import pytest
import bpy
from .constants import data_dir

cellpack_dir = data_dir / "cellpack"


files_to_test = [f for f in cellpack_dir.glob("*")]


@pytest.mark.parametrize("file", files_to_test)
def test_load_petworld(file):
    ens = mn.entities.ensemble.load_cellpack(
        file_path=file,
        name="CellPack",
        node_setup=False,
        fraction=0.1,
    )
    assert False


# @pytest.mark.parametrize("format", ["bcif", "cif"])
# def test_load_cellpack(snapshot, format):
#     name = f"Cellpack_{format}"
#     ens = mn.entities.ensemble.load_cellpack(
#         data_dir / f"square1.{format}", name=name, node_setup=False, fraction=0.1
#     )

#     assert ens.name == name
#     assert snapshot == str(ens.object["chain_ids"])
#     obj_names = [obj.name for obj in ens.instance_collection.objects]
#     assert snapshot == "\n".join(obj_names)

#     ens.node_group.nodes["Ensemble Instance"].inputs["As Points"].default_value = False
#     mn.blender.nodes.realize_instances(ens.object)
#     for attribute in ens.list_attributes():
#         assert snapshot == ens.named_attribute(attribute)

#     pos_eval = ens.named_attribute("position", evaluate=True)
#     assert snapshot == pos_eval.shape
#     assert snapshot == pos_eval

#     obj = ens.object
