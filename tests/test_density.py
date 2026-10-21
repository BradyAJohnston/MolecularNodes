import itertools
import bpy
import databpy
import numpy as np
import pytest
from databpy import ObjectTracker
import molecularnodes as mn
from molecularnodes.nodes._utils import get_final_style_nodes
from .constants import data_dir
from .utils import NumpySnapshotExtension


@pytest.fixture
def density_file(isolated_density_file):
    file = isolated_density_file(data_dir / "emd_24805.map.gz")
    # Make all densities are removed
    for o in bpy.data.objects:
        if o.mn.entity_type == "density":
            bpy.data.objects.remove(o, do_unlink=True)
    return file


def test_density_load(density_file):
    density = mn.entities.density.Grids.load(density_file)
    pos = density.named_attribute("position")

    assert len(pos) > 1000

    avg = np.mean(pos, axis=0)
    assert np.linalg.norm(avg) > 0.5

    print(f"{list(bpy.data.objects)=}")

    assert density.props.entity_type == "density"
    assert density.object.users_collection[0] == mn.blender.coll.mn()


def test_density_centered(density_file):
    # First load using standard parameters to test recreation of vdb
    density = mn.entities.density.Grids.load(density_file, center=True, overwrite=True)
    pos = density.named_attribute("position")
    avg = np.mean(pos, axis=0)
    assert len(pos) > 1000
    # centered on the origin, vs ~20.0 for the same volume uncentered
    assert np.linalg.norm(avg) < 1.0


def test_density_invert(density_file):
    # First load using standard parameters to test recreation of vdb
    density = mn.entities.density.Grids.load(density_file)

    # Then refresh the scene
    bpy.data.objects.remove(density.object, do_unlink=True)

    density = mn.entities.density.Grids.load(density_file, invert=True)
    style_node = get_final_style_nodes(density.modifier_node_tree)[0]
    style_node.inputs["Threshold"].default_value = 0.01

    pos = density.named_attribute("position")
    # At this threshold after inverting we should have a cube the size of the volume
    assert pos[:, 0].max() > 2.0
    assert pos[:, 1].max() > 2.0
    assert pos[:, 2].max() > 2.0


def test_density_multiple_load(isolated_density_file):
    file = isolated_density_file(data_dir / "emd_24805.map.gz")
    density1 = mn.entities.density.Grids.load(file)
    density2 = mn.entities.density.Grids.load(file)

    assert density1.props.entity_type == "density"
    assert density2.props.entity_type == "density"
    assert density1.object.users_collection[0] == mn.blender.coll.mn()
    assert density2.object.users_collection[0] == mn.blender.coll.mn()


def test_density_naming_op(density_file):
    bpy.ops.mn.import_density(filepath=str(density_file))
    _object = bpy.data.objects[density_file.name]


@pytest.mark.parametrize(
    "invert,center", list(itertools.product([True, False], repeat=2))
)
def test_density_operator(
    snapshot_custom: NumpySnapshotExtension, density_file, invert, center
):
    scene = bpy.context.scene
    with ObjectTracker() as o:
        bpy.ops.mn.import_density(
            filepath=str(density_file),
            invert=invert,
            center=center,
        )
        density: mn.entities.Density = scene.MNSession.match(o.latest())

    n_vertices = len(databpy.GeometrySet(density.object).named_attribute("position"))
    assert n_vertices == [198862, 287422][int(invert)]


@pytest.fixture
def density_file_dx(isolated_density_file):
    file = isolated_density_file(data_dir / "water.dx.gz")
    # Make sure all densities are removed
    for o in bpy.data.objects:
        if o.mn.entity_type == "density":
            bpy.data.objects.remove(o, do_unlink=True)
    return file


def test_density_load_dx(density_file_dx):
    density = mn.entities.density.Grids.load(density_file_dx)
    print(f"{list(bpy.data.objects)=}")
    assert density.props.entity_type == "density"
    assert density.object.users_collection[0] == mn.blender.coll.mn()


# this test fails without the fallback using mrcfile
def test_fallback_reading(isolated_density_file):
    mn.entities.density.Grids.load(
        isolated_density_file(data_dir / "62270-small_sg0.mrc")
    )


def test_fallback_transforms(isolated_density_file):
    def evaluated_positions(file):
        density = mn.entities.density.Grids.load(file, overwrite=True)
        return databpy.GeometrySet(density.object).named_attribute("position")

    pos_default = evaluated_positions(
        isolated_density_file(data_dir / "62270-small_sg1.mrc")
    )
    pos_chimerax = evaluated_positions(
        isolated_density_file(data_dir / "62270-small_sg0.mrc")
    )

    assert np.allclose(pos_default, pos_chimerax)


# signed grids, affine transforms and sampling a grid ------------------------------


def _write_dx(path, data, origin, axes):
    """Write a .dx file with a full `delta` matrix (rows are the voxel axes, in Å)."""
    nx, ny, nz = data.shape
    values = np.asarray(data, dtype=float).ravel()  # z fastest, as OpenDX expects
    lines = [
        f"object 1 class gridpositions counts {nx} {ny} {nz}",
        "origin " + " ".join(f"{v:.8g}" for v in origin),
        *("delta " + " ".join(f"{v:.8g}" for v in row) for row in axes),
        f"object 2 class gridconnections counts {nx} {ny} {nz}",
        f"object 3 class array type double rank 0 items {values.size} data follows",
        *(
            " ".join(f"{v:.10g}" for v in values[i : i + 3])
            for i in range(0, values.size, 3)
        ),
        'object "density" class field',
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


def _linear_field(shape):
    i, j, k = np.indices(shape, dtype=float)
    # trilinear interpolation reproduces a linear field exactly
    return i + 10 * j + 100 * k


def _input_default(group_name, socket_name):
    tree = bpy.data.node_groups[group_name]
    return np.array(
        next(
            item.default_value
            for item in tree.interface.items_tree
            if item.item_type == "SOCKET"
            and item.in_out == "INPUT"
            and item.name == socket_name
        )
    )


def _sample_grid_at(density, points_angstrom):
    """Evaluate `Color Sample Grid` at points (Å) on a separate mesh object."""
    from nodebpy import TreeBuilder
    from nodebpy import geometry as g
    from molecularnodes.nodes.geometry import ColorSampleGrid

    points = np.asarray(points_angstrom, dtype=float) * 0.1
    obj = bpy.data.objects.new("sample_points", bpy.data.meshes.new("sample_points"))
    bpy.context.scene.collection.objects.link(obj)
    obj.data.from_pydata(points.tolist(), [], [])
    tree = bpy.data.node_groups.new("sample_points", "GeometryNodeTree")
    obj.modifiers.new("sample", "NODES").node_group = tree
    with TreeBuilder(tree) as t:
        geometry = t.inputs.geometry("Geometry")
        result = t.outputs.geometry("Geometry")
        node = ColorSampleGrid(object=density.object)
        (
            geometry
            >> g.StoreNamedAttribute.point.float(name="value", value=node.o.value)
            >> g.StoreNamedAttribute.point.color(name="sampled", value=node.o.color)
            >> result
        )
    geometry_set = databpy.GeometrySet(obj)
    return (
        geometry_set.named_attribute("value"),
        geometry_set.named_attribute("sampled"),
    )


@pytest.mark.parametrize("threshold", [1.5, -1.5])
def test_iso_surface_signed_threshold(tmp_path, threshold):
    from nodebpy import geometry as g
    from molecularnodes.nodes.geometry import DensityStyleISOSurface

    # +3 and -3 Gaussian lobes centred at x = +4 and x = -4 Å
    x = (np.arange(41) - 20) * 0.5
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    centres = {"positive": np.array([4.0, 0, 0]), "negative": np.array([-4.0, 0, 0])}

    def blob(c):
        return np.exp(-((X - c[0]) ** 2 + (Y - c[1]) ** 2 + (Z - c[2]) ** 2) / 4.0)

    data = 3.0 * blob(centres["positive"]) - 3.0 * blob(centres["negative"])
    path = _write_dx(tmp_path / "signed.dx", data, [x[0]] * 3, np.eye(3) * 0.5)
    density = mn.entities.density.Grids.load(path, style=None)
    with density.tree.reset(input="Volume") as (volume, join):
        (
            volume
            >> DensityStyleISOSurface(threshold=threshold)
            >> g.StoreNamedAttribute.point.vector(name="normal", value=g.Normal())
            >> join
        )
    # outside the old 0-1 clamp
    style = next(
        node
        for node in density.modifier_node_tree.nodes
        if node.bl_idname == "GeometryNodeGroup"
        and node.node_tree.name.startswith(DensityStyleISOSurface._name)
    )
    assert style.inputs["Threshold"].default_value == pytest.approx(threshold)

    positions = density.named_attribute("position") / 0.1
    normals = density.named_attribute("normal")
    # 3 exp(-r^2 / 4) = 1.5 at r = sqrt(4 ln 2) = 1.67 Å; a threshold clamped to 1
    # would give sqrt(4 ln 3) = 2.10 Å
    expected_radius = np.sqrt(4 * np.log(2))
    for side, selection in [
        ("positive", positions[:, 0] > 0),
        ("negative", positions[:, 0] < 0),
    ]:
        assert selection.sum() > 50, side
        radial = positions[selection] - centres[side]
        distance = np.linalg.norm(radial, axis=1)
        assert distance.mean() == pytest.approx(expected_radius, abs=0.05), side
        # both lobes are closed surfaces with outward normals
        outward = np.sum(radial / distance[:, None] * normals[selection], axis=1)
        assert outward.min() > 0.9, side


def test_density_orthogonal_transform_unchanged(tmp_path):
    shape = (6, 7, 8)
    origin = np.array([1.5, -2.0, 3.25])
    delta = np.array([0.5, 0.75, 1.25])
    path = _write_dx(
        tmp_path / "ortho.dx", _linear_field(shape), origin, np.diag(delta)
    )
    density = mn.entities.density.Grids.load(path, style=None)
    assert np.array_equal(density.axes, np.diag(delta))

    density.object.data.grids.load()
    grid = density.object.data.grids["density"]
    expected = np.eye(4)
    expected[:3, :3] = np.diag(delta * 0.1)
    expected[:3, 3] = origin * 0.1
    assert np.allclose(np.array(grid.matrix_object), expected, atol=1e-7)

    bpy.utils.expose_bundled_modules()
    import openvdb as vdb  # type: ignore

    # written with scale + translate as before, not the affine path
    metadata = vdb.readAllGridMetadata(density.file_vdb)[0]
    assert "MN_affine" not in metadata
    del metadata

    ijk = np.array([[1.5, 2.25, 3.0], [4.0, 4.0, 4.0], [0.0, 6.0, 2.5]])
    values, _ = _sample_grid_at(density, origin + ijk * delta)
    assert np.allclose(values, ijk @ [1.0, 10.0, 100.0], atol=1e-3)


def test_density_affine_transform_dx(tmp_path):
    shape = (6, 7, 8)
    origin = np.array([2.0, -1.0, 3.0])
    # sheared voxel axes, as rows
    axes = np.array([[1.0, 0.0, 0.0], [0.5, 1.0, 0.0], [0.2, 0.3, 1.2]])
    path = _write_dx(tmp_path / "sheared.dx", _linear_field(shape), origin, axes)
    density = mn.entities.density.Grids.load(path, style=None)
    assert np.allclose(density.axes, axes)

    density.object.data.grids.load()
    grid = density.object.data.grids["density"]
    expected = np.eye(4)
    expected[:3, :3] = axes.T * 0.1  # Blender's matrix holds the axes as columns
    expected[:3, 3] = origin * 0.1
    assert np.allclose(np.array(grid.matrix_object), expected, atol=1e-7)

    # the world position of a (fractional) voxel index holds the value at that index
    ijk = np.array([[1.5, 2.25, 3.0], [4.0, 4.0, 4.0], [0.5, 6.0, 2.0], [0, 0, 0]])
    values, _ = _sample_grid_at(density, origin + ijk @ axes)
    assert np.allclose(values, ijk @ [1.0, 10.0, 100.0], atol=1e-3)


def test_density_affine_transform_mrc(isolated_density_file):
    # monoclinic map (beta = 94.3 deg) with axis order 3, 1, 2 and non-zero start
    # indices, shipped with the GridDataFormats tests
    from pathlib import Path
    import mrcfile
    from gridData.tests import datafiles

    path = isolated_density_file(Path(datafiles.__file__).parent / "EMD-3001.map.bz2")
    with mrcfile.open(path, permissive=True) as mrc:
        h = mrc.header
        data = mrc.data.copy()  # indexed [section, row, column]
        lengths = np.array([h.cella.x, h.cella.y, h.cella.z], dtype=float)
        alpha, beta, gamma = np.radians([h.cellb.alpha, h.cellb.beta, h.cellb.gamma])
        sampling = np.array([h.mx, h.my, h.mz], dtype=float)
        # crystal axis along columns, rows and sections, and their start indices
        order = np.array([h.mapc, h.mapr, h.maps]) - 1
        starts = np.array([h.nxstart, h.nystart, h.nzstart], dtype=float)
    assert not np.isclose(np.degrees(beta), 90.0)

    # fractional -> Cartesian with a along x and b in the xy plane
    cy = (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cell = np.array(
        [
            [lengths[0], 0, 0],
            [lengths[1] * np.cos(gamma), lengths[1] * np.sin(gamma), 0],
            [
                lengths[2] * np.cos(beta),
                lengths[2] * cy,
                lengths[2] * np.sqrt(1 - np.cos(beta) ** 2 - cy**2),
            ],
        ]
    )

    density = mn.entities.density.Grids.load(path, style=None)

    rng = np.random.default_rng(0)
    indices = np.column_stack([rng.integers(0, n, 12) for n in data.shape])
    points = []
    for s, r, c in indices:
        grid_index = np.zeros(3)
        grid_index[order] = np.array([c, r, s]) + starts
        points.append((grid_index / sampling) @ cell)
    values, _ = _sample_grid_at(density, points)
    expected = data[indices[:, 0], indices[:, 1], indices[:, 2]]
    assert np.allclose(values, expected, atol=1e-5)


def test_color_sample_grid_values_and_colors(tmp_path):
    from molecularnodes.nodes.geometry import ColorSampleGrid

    # value -20 + 10 i along x, on a 1 Å grid at the origin
    data = np.broadcast_to((np.arange(5) * 10.0 - 20.0)[:, None, None], (5, 3, 3))
    path = _write_dx(tmp_path / "ramp.dx", data, [0, 0, 0], np.eye(3))
    density = mn.entities.density.Grids.load(path, style=None)
    points = [[0, 1, 1], [1.75, 1, 1], [2, 1, 1], [2.25, 1, 1], [4, 1, 1], [20, 1, 1]]
    values, colors = _sample_grid_at(density, points)
    # trilinear inside the grid, 0 outside it
    assert np.allclose(values, [-20, -2.5, 0, 2.5, 20, 0], atol=1e-4)

    name = ColorSampleGrid._name
    assert _input_default(name, "Min") == -5.0
    assert _input_default(name, "Max") == 5.0
    a, mid, b = (_input_default(name, n)[:3] for n in ("A", "Intermediate", "B"))
    # clamped below Min and above Max, linear through the Intermediate color
    expected = [a, (a + mid) / 2, mid, (mid + b) / 2, b, mid]
    assert np.allclose(colors[:, :3], expected, atol=1e-4)


def test_add_style_color_by_density(fetch, tmp_path):
    from molecularnodes.nodes.geometry import ColorSampleGrid

    mol = fetch("1cbs")
    positions = mol.named_attribute("position") / mol.world_scale
    origin = positions.min(axis=0) - 10.0
    shape = tuple(np.ceil(np.ptp(positions, axis=0) + 20).astype(int))
    path = _write_dx(tmp_path / "positive.dx", np.full(shape, 10.0), origin, np.eye(3))
    density = mn.entities.density.Grids.load(path, style=None)

    mol.add_style("surface", color=density)
    colors = mol.named_attribute("Color", evaluate=True)
    # every surface vertex is in the +10 grid, above Max, so gets color B
    assert len(colors) > 1000
    b = _input_default(ColorSampleGrid._name, "B")
    assert np.allclose(colors[:, :3], b[:3], atol=1e-4)


def _ramp_around(mol, tmp_path, name):
    "A 1 Å grid around the molecule whose value is x (Å) minus the molecule's mean x."
    positions = mol.named_attribute("position") / mol.world_scale
    origin = np.floor(positions.min(axis=0)) - 10.0
    shape = tuple(np.ceil(np.ptp(positions, axis=0) + 20).astype(int))
    x0 = float(positions[:, 0].mean())
    data = np.broadcast_to(
        (origin[0] + np.arange(shape[0]) - x0)[:, None, None], shape
    ).copy()
    return _write_dx(tmp_path / name, data, origin, np.eye(3)), x0


def test_density_grid_object(tmp_path):
    path = _write_dx(tmp_path / "grid.dx", np.ones((4, 4, 4)), [0, 0, 0], np.eye(3))
    density = mn.entities.density.Grids.load(path, style="density_surface")
    grid_object = density.grid_object

    assert grid_object.name.startswith(".data_")
    assert grid_object.data == density.object.data
    assert grid_object.parent == density.object
    assert len(grid_object.modifiers) == 0
    assert grid_object.hide_render and grid_object.hide_viewport
    assert grid_object.users_collection[0] == mn.blender.coll.data()
    # not drawn or rendered: its collection is excluded from the view layer
    assert grid_object not in list(bpy.context.view_layer.objects)
    # the same object on every access
    assert density.grid_object == grid_object
    # while the styled density itself no longer holds a grid
    assert databpy.GeometrySet(density.object).volume is None


def test_density_grid_object_removed_with_density(tmp_path):
    path = _write_dx(tmp_path / "grid.dx", np.ones((4, 4, 4)), [0, 0, 0], np.eye(3))
    density = mn.entities.density.Grids.load(path)
    name = density.grid_object.name
    mn.session.get_session().remove(density.uuid)
    assert name not in bpy.data.objects


def test_add_style_color_by_styled_density(fetch, tmp_path):
    from nodebpy import geometry as g
    from molecularnodes.nodes.geometry import ColorSampleGrid, StyleSurface

    mol = fetch("1cbs")
    path, x0 = _ramp_around(mol, tmp_path, "ramp.dx")
    # styled, so its own object holds a mesh rather than a grid
    density = mn.entities.density.Grids.load(path, style="density_iso_surface")

    with mol.tree.reset() as (atoms, join):
        (
            atoms
            >> StyleSurface()
            >> g.StoreNamedAttribute.point.float(
                name="value",
                value=ColorSampleGrid(object=density.grid_object).o.value,
            )
            >> join
        )
    surface = mol.named_attribute("position", evaluate=True) / mol.world_scale
    values = mol.named_attribute("value", evaluate=True)
    assert len(values) > 1000
    assert np.allclose(values, surface[:, 0] - x0, atol=1e-3)

    with mol.tree.reset():
        pass
    mol.add_style("surface", color=density)
    surface = mol.named_attribute("position", evaluate=True) / mol.world_scale
    colors = mol.named_attribute("Color", evaluate=True)[:, :3]
    a, mid, b = (
        _input_default(ColorSampleGrid._name, n)[:3] for n in "A Intermediate B".split()
    )
    factor = np.clip((surface[:, 0] - x0 + 5.0) / 10.0, 0, 1)[:, None]
    low = a + np.clip(2 * factor, 0, 1) * (mid - a)
    expected = low + np.clip(2 * factor - 1, 0, 1) * (b - low)
    assert np.allclose(colors, expected, atol=1e-3)
    # the ramp spans the -5 to 5 range, so both ends are used
    assert (factor == 0).any() and (factor == 1).any()
