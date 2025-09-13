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
