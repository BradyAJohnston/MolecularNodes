import gzip
import shutil
from pathlib import Path
import pytest
import molecularnodes as mn
from molecularnodes.nodes import nodes
from .constants import data_dir

cellpack_dir = data_dir / "cellpack/petworld"


files_to_test = [f for f in cellpack_dir.glob("*") if str(f).endswith((".gz", ".bcif"))]


def maybe_unzip(file):
    if str(file).endswith(".gz"):
        # Create a temporary unzipped file
        unzipped_path = str(file)[:-3]  # Remove .gz extension
        if not Path(unzipped_path).exists():
            with gzip.open(file, "rb") as f_in:
                with open(unzipped_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
        return unzipped_path
    else:
        return file


@pytest.mark.parametrize("file", files_to_test)
def test_load_petworld(file):
    file_path = maybe_unzip(data_dir / "cellpack" / file)

    _ens = mn.entities.ensemble.load_cellpack(
        file_path=file_path,
        name="CellPack",
        node_setup=False,
        fraction=0.1,
    )


@pytest.mark.parametrize("format", ["bcif", "cif"])
def test_load_cellpack(snapshot, format):
    file_path = data_dir / f"cellpack/square1.{format}"

    ens = mn.entities.ensemble.load_cellpack(file_path, node_setup=False, fraction=0.1)

    assert ens.name == Path(file_path).name
    assert snapshot == str(ens.object["chain_ids"])
    obj_names = [obj.name for obj in ens.instance_collection.objects]
    assert snapshot == "\n".join(obj_names)

    ens.node_group.nodes["Ensemble Instance"].inputs["As Points"].default_value = False
    nodes.realize_instances(ens.object)
    for attribute in ens.list_attributes():
        assert snapshot == ens.named_attribute(attribute)

    pos_eval = ens.named_attribute("position", evaluate=True)
    assert snapshot == pos_eval.shape
    assert snapshot == pos_eval

    _obj = ens.object
