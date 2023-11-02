import molecularnodes as mn

def test_coll():
    assert mn.coll.mn().name == "MolecularNodes"
    assert mn.coll.mn().name == "MolecularNodes"
    assert mn.coll.data().name == "MN_data"
    assert mn.coll.cellpack().name == "cellpack_"
    assert mn.coll.cellpack(fallback=True).name == "cellpack_"
    assert mn.coll.cellpack().name == "cellpack_.001"
    assert mn.coll.frames().name == "_frames"
    assert mn.coll.frames("4OZS").name == "4OZS_frames"