import sys
import types
import numpy as np


def _make_test_mrc(
    path, shape=(4, 4, 4), voxel_size=(1.5, 2.0, 2.5), origin=(10.0, 20.0, 30.0)
):
    import mrcfile

    data = np.arange(np.prod(shape), dtype=np.float32).reshape(shape)
    with mrcfile.new(path, overwrite=True) as mrc:
        mrc.set_data(data)
        # Set voxel size if available
        try:
            mrc.voxel_size = voxel_size
        except Exception:
            # fall back to header cell dimensions
            nx, ny, nz = shape[2], shape[1], shape[0]
            ax, ay, az = voxel_size[0] * nx, voxel_size[1] * ny, voxel_size[2] * nz
            mrc.header.cella.x = ax
            mrc.header.cella.y = ay
            mrc.header.cella.z = az
        # Set origin if supported
        try:
            mrc.header.origin.x = origin[0]
            mrc.header.origin.y = origin[1]
            mrc.header.origin.z = origin[2]
        except Exception:
            pass
    return (
        np.asarray(data),
        np.asarray(voxel_size, dtype=float),
        np.asarray(origin, dtype=float),
    )


def test_parse_grid_with_mrcfile_fallback(monkeypatch, tmp_path):
    # Provide lightweight stubs for modules imported by grids.py
    if "bpy" not in sys.modules:
        bpy_stub = types.SimpleNamespace()
        bpy_stub.utils = types.SimpleNamespace(expose_bundled_modules=lambda: None)
        bpy_stub.app = types.SimpleNamespace(version=(4, 0, 0))
        sys.modules["bpy"] = bpy_stub
    if "databpy" not in sys.modules:
        sys.modules["databpy"] = types.SimpleNamespace()
    if "gridData" not in sys.modules:
        griddata_stub = types.ModuleType("gridData")

        class _DummyGrid:  # placeholder; we will monkeypatch to raising below
            def __init__(self, *args, **kwargs):
                pass

        griddata_stub.Grid = _DummyGrid
        sys.modules["gridData"] = griddata_stub

    from molecularnodes.entities.density import grids as grids_mod

    # Create a sample .mrc file
    mrc_path = tmp_path / "sample.mrc"
    data, voxel_size, origin = _make_test_mrc(str(mrc_path))

    # Force GridDataFormats.Grid to raise to exercise the fallback
    class _RaisingGrid:
        def __init__(self, *args, **kwargs):
            raise RuntimeError("Simulated GridDataFormats parse failure")

    monkeypatch.setattr(grids_mod, "Grid", _RaisingGrid)

    # Call the fallback parser directly to avoid openvdb/blender dependencies
    g = object.__new__(grids_mod.Grids)  # bypass __init__ which writes VDB
    md = {"filepath": str(mrc_path), "invert": False, "center": False}
    parsed = grids_mod.Grids._parse_grid_with_fallback(g, str(mrc_path), None, md)

    # Validate parsed content
    assert (
        hasattr(parsed, "grid")
        and hasattr(parsed, "delta")
        and hasattr(parsed, "origin")
    )
    assert parsed.grid.shape == data.shape
    assert np.allclose(parsed.delta, voxel_size)
    assert np.allclose(parsed.origin, origin)

    # Ensure metadata is carried through
    assert getattr(parsed, "metadata", {}) == md
