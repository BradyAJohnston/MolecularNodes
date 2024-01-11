import molecularnodes as mn


def test_coll():
    assert mn.blender.coll.mn().name == "MolecularNodes"
    assert mn.blender.coll.mn().name == "MolecularNodes"
    assert mn.blender.coll.data().name == "MN_data"
    assert mn.blender.coll.cellpack().name == "cellpack_"
    assert mn.blender.coll.cellpack(fallback=True).name == "cellpack_"
    assert mn.blender.coll.cellpack().name == "cellpack_.001"
    assert mn.blender.coll.frames().name == "_frames"
    assert mn.blender.coll.frames("4OZS").name == "4OZS_frames"
