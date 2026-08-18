import gzip
import shutil
from pathlib import Path
import pytest
import molecularnodes as mn
from molecularnodes.nodes import nodes
from .constants import data_dir
from .utils import round_significant

cellpack_dir = data_dir / "cellpack/petworld"


files_to_test = [f for f in cellpack_dir.glob("*") if f.suffix in (".gz", ".bcif")]


def maybe_unzip(file: Path) -> Path:
    if file.suffix == ".gz":
        unzipped_path = file.with_suffix("")
        if not unzipped_path.exists():
            with gzip.open(file, "rb") as f_in:
                with open(unzipped_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
        return unzipped_path
    else:
        return file


@pytest.mark.parametrize("file", files_to_test)
def test_load_petworld(file):
    file_path = maybe_unzip(data_dir / "cellpack" / file)

    _ens = mn.entities.ensemble.CellPack.load(
        file_path=file_path,
        name="CellPack",
        node_setup=False,
    )


@pytest.mark.parametrize("format", ["bcif", "cif"])
def test_load_cellpack(snapshot, format):
    file_path = data_dir / f"cellpack/square1.{format}"

    ens = mn.entities.ensemble.CellPack.load(file_path, node_setup=False)
    assert ens._entity_type == mn.entities.base.EntityType.ENSEMBLE_CELLPACK
    assert ens.props.entity_type == ens._entity_type.value

    assert ens.name == Path(file_path).name
    assert snapshot == str(ens.props.chain_ids)
    obj_names = [obj.name for obj in ens.instance_collection.objects]
    assert snapshot == "\n".join(obj_names)

    ens.node_group.nodes["Ensemble Instance"].inputs["As Points"].default_value = False
    nodes.realize_instances(ens.object)
    for attribute in ens.list_attributes():
        assert snapshot == ens[attribute]

    pos_eval = ens.named_attribute("position", evaluate=True)
    assert snapshot == pos_eval.shape
    # positions evaluated through geometry nodes carry platform-specific last-bit
    # differences, so quantise before comparing
    assert snapshot == round_significant(pos_eval)

    _obj = ens.object
