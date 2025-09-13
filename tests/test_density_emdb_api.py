import numpy as np
import pytest
import molecularnodes as mn
from .emdb_pooch import fetch_emdb_map


@pytest.fixture
def emdb_density_map():
    file = fetch_emdb_map("EMD-48397")
    # we don't just use Path().stem because usually is .map.gz
    base = file.name.split(".")[0]
    vdb_path = file.parent / f"{base}.vdb"
    vdb_path.unlink(missing_ok=True)
    return file


def test_emdb_api_density_load(emdb_density_map):
    density = mn.entities.density.load(emdb_density_map)
    pos = density.named_attribute("position")

    # Basic sanity checks only; remote data can change slightly over time
    assert len(pos) > 1000
    avg = np.mean(pos, axis=0)
    assert np.linalg.norm(avg) > 0.1

    assert density.object.mn.entity_type == "density"
    assert density.object.users_collection[0] == mn.blender.coll.mn()
