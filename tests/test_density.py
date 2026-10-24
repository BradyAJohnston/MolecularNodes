import itertools
from pathlib import Path
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


def _grid_matrix(density, name="density"):
    """The 4x4 transform Blender gives a density's grid, read with `Grid Info`."""
    from nodebpy import TreeBuilder
    from nodebpy import geometry as g

    obj = bpy.data.objects.new("grid_matrix", bpy.data.meshes.new("grid_matrix"))
    bpy.context.scene.collection.objects.link(obj)
    tree = bpy.data.node_groups.new("grid_matrix", "GeometryNodeTree")
    obj.modifiers.new("probe", "NODES").node_group = tree
    with TreeBuilder(tree) as t:
        result = t.outputs.geometry("Geometry")
        grid = g.GetBundleItem.float(
            g.ObjectInfo(density.object).o.geometry.bundle(),
            f"MN Grid {name}",
            structure_type="GRID",
        )
        (
            g.Points(1)
            >> g.StoreNamedAttribute.point.matrix(
                name="matrix", value=g.GridInfo(grid.o.item).o.transform
            )
            >> result
        )
    return databpy.GeometrySet(obj).named_attribute("matrix")[0]


def _sample_grid_at(density, points_angstrom, grid_name="density"):
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
        node = ColorSampleGrid(object=density.object, grid_name=grid_name)
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
    with density.tree.reset() as (_, join):
        (
            density.import_volume()
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

    expected = np.eye(4)
    expected[:3, :3] = np.diag(delta * 0.1)
    expected[:3, 3] = origin * 0.1
    assert np.allclose(_grid_matrix(density), expected, atol=1e-7)

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

    expected = np.eye(4)
    expected[:3, :3] = axes.T * 0.1  # Blender's matrix holds the axes as columns
    expected[:3, 3] = origin * 0.1
    assert np.allclose(_grid_matrix(density), expected, atol=1e-7)

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


def _ramp_density(tmp_path, style, shape=(8, 8, 8)):
    "A linear field on a 1 Å grid at the origin, loaded with the given style."
    path = _write_dx(tmp_path / "ramp.dx", _linear_field(shape), [0, 0, 0], np.eye(3))
    return mn.entities.density.Grids.load(path, style=style)


RAMP_POINTS = np.array([[1.5, 2.25, 3.0], [4.0, 4.0, 4.0], [0.5, 6.0, 2.0]])
RAMP_VALUES = RAMP_POINTS @ [1.0, 10.0, 100.0]


def test_styled_density_keeps_volume_in_bundle(tmp_path):
    density = _ramp_density(tmp_path, "density_iso_surface")
    # the styled output is a mesh: the bundled volume is not a rendered component
    geometry = databpy.GeometrySet(density.object)
    assert geometry.volume is None
    assert geometry.mesh is not None
    # but the grid is still there to sample, through the `MN Grid density` bundle item
    values, _ = _sample_grid_at(density, RAMP_POINTS)
    assert np.allclose(values, RAMP_VALUES, atol=1e-3)


def test_unstyled_density_samples_its_own_volume(tmp_path):
    density = _ramp_density(tmp_path, None)
    assert databpy.GeometrySet(density.object).volume is not None
    values, _ = _sample_grid_at(density, RAMP_POINTS)
    assert np.allclose(values, RAMP_VALUES, atol=1e-3)


def test_density_styles_joined_keep_volume(tmp_path):
    from molecularnodes.nodes.geometry import (
        DensityStyleISOSurface,
        DensityStyleSurface,
        ImportVolume,
    )

    density = _ramp_density(tmp_path, None)
    with density.tree.reset() as (_, join):
        density.import_volume() >> DensityStyleISOSurface(threshold=300.0) >> join
        density.import_volume() >> DensityStyleSurface(threshold=500.0) >> join
    # a third style added to the existing tree goes into the same join
    with density.tree as tree:
        density.import_volume() >> DensityStyleSurface(threshold=200.0) >> tree.join

    nodes = density.modifier_node_tree.nodes
    imports = [
        n
        for n in nodes
        if n.bl_idname == "GeometryNodeGroup"
        and n.node_tree.name.startswith(ImportVolume._name)
    ]
    joins = [n for n in nodes if n.bl_idname == "GeometryNodeJoinGeometry"]
    assert len(imports) == 3 and len(joins) == 1
    assert len(joins[0].inputs[0].links) == 3
    assert len(get_final_style_nodes(density.modifier_node_tree)) == 3

    assert databpy.GeometrySet(density.object).volume is None
    values, _ = _sample_grid_at(density, RAMP_POINTS)
    assert np.allclose(values, RAMP_VALUES, atol=1e-3)


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
                value=ColorSampleGrid(object=density.object).o.value,
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


def test_apbs_1brs_example_grid(isolated_density_file):
    "The barnase-barstar potential used in the electrostatics docs page."
    from gridData import Grid

    path = isolated_density_file(data_dir / "1brs_AD_apbs.dx.gz")
    density = mn.entities.density.Grids.load(path, style="density_iso_surface")
    assert density.grid.grid.shape == (65, 65, 65)

    # side-chain charge centres from PDB 1BRS (Å): barnase (chain A) lysines and
    # barstar (chain D) carboxylates
    atoms = {
        "A LYS27 NZ": (28.301, 44.689, 7.816),
        "A LYS62 NZ": (22.219, 30.830, 2.640),
        "A LYS108 NZ": (6.643, 38.397, 5.241),
        "D ASP39 OD1": (28.063, 40.822, 10.213),
        "D GLU76 OE1": (29.816, 44.208, -2.210),
    }
    points = np.array(list(atoms.values()))
    values, _ = _sample_grid_at(density, points)

    # trilinear interpolation in numpy, as Sample Grid does
    reference = Grid(str(path))
    index = (points - reference.origin) / reference.delta
    base = np.floor(index).astype(int)
    frac = index - base
    expected = np.zeros(len(points))
    for corner in np.ndindex(2, 2, 2):
        weight = np.prod(np.where(corner, frac, 1 - frac), axis=1)
        i, j, k = (base + corner).T
        expected += weight * reference.grid[i, j, k]
    # values change by ~50 kT/e per Å at these atoms, so allow for float32 positions
    assert np.allclose(values, expected, atol=1.0)
    for name, value in zip(atoms, values):
        assert (value > 10) if "LYS" in name else (value < -10), (name, value)


# importing volumes in the node tree --------------------------------------------------


def _constant_dx(path, value, shape=(4, 4, 4)):
    return _write_dx(path, np.full(shape, float(value)), [0, 0, 0], np.eye(3))


def test_density_object_is_empty_mesh(tmp_path):
    # the volume is imported by the node tree, not by an operator (#893)
    density = _ramp_density(tmp_path, "density_surface")
    assert density.object.type == "MESH"
    assert len(density.object.data.vertices) == 0
    assert len(density.named_attribute("position")) > 0


def test_density_tree_has_no_input(tmp_path):
    density = _ramp_density(tmp_path, None)
    with density.tree.reset() as (unused, join):
        assert unused is None
        density.import_volume() >> join
    inputs = [
        item
        for item in density.modifier_node_tree.interface.items_tree
        if item.item_type == "SOCKET" and item.in_out == "INPUT"
    ]
    assert inputs == []
    with pytest.raises(AttributeError, match="ImportVolume"):
        density.tree.atoms


def test_density_sequence_follows_frame(tmp_path):
    # several files load into one density and play one per frame (#1100)
    paths = [_constant_dx(tmp_path / f"frame_{i}.dx", 10 * (i + 1)) for i in range(3)]
    density = mn.entities.density.Grids.load(paths, style=None, frame_start=5)
    assert len(density.files_vdb) == 3
    assert all(Path(f).exists() for f in density.files_vdb)
    assert density.vdb_path.endswith("_{index}.vdb")

    point = [[1.5, 1.5, 1.5]]
    # before the start holds the first file, past the end holds the last
    for frame, expected in [(1, 10), (5, 10), (6, 20), (7, 30), (20, 30)]:
        bpy.context.scene.frame_set(frame)
        values, _ = _sample_grid_at(density, point)
        assert values[0] == pytest.approx(expected), frame


def test_density_sequence_by_index(tmp_path):
    paths = [_constant_dx(tmp_path / f"grid_{i}.dx", i + 1) for i in range(2)]
    density = mn.entities.density.Grids.load(paths, style=None, follow_frame=False)
    # two branches import different files, kept apart by their names
    with density.tree.reset() as (_, join):
        density.import_volume(name="first", index=0) >> join
        density.import_volume(name="second", index=1) >> join
    point = [[1.5, 1.5, 1.5]]
    for name, expected in [("first", 1.0), ("second", 2.0)]:
        values, _ = _sample_grid_at(density, point, grid_name=name)
        assert values[0] == pytest.approx(expected), name


def test_density_path_with_braces(tmp_path):
    folder = tmp_path / "{index} dir"
    folder.mkdir()
    density = _ramp_density(folder, None)
    assert density.vdb_path == density.file_vdb.replace("{", "{{").replace("}", "}}")
    values, _ = _sample_grid_at(density, RAMP_POINTS)
    assert np.allclose(values, RAMP_VALUES, atol=1e-3)


# Gaussian cube files ---------------------------------------------------------------


def _write_cube(path, data, origin, axes, angstrom=False, dataset_ids=None):
    """Write a cube file; `data` is (n_datasets, n1, n2, n3), lengths in file units."""
    data = np.asarray(data, dtype=float)
    shape = data.shape[1:]
    sign = -1 if angstrom else 1
    n_atoms = -2 if dataset_ids else 2
    lines = [
        "comment one",
        "comment two",
        f"{n_atoms} " + " ".join(f"{v:.8f}" for v in origin),
        *(
            f"{sign * n} " + " ".join(f"{v:.8f}" for v in row)
            for n, row in zip(shape, axes)
        ),
        "8 8.0 0.0 0.0 0.0",
        "1 1.0 1.5 0.0 0.0",
    ]
    if dataset_ids:
        lines.append(" ".join(str(x) for x in [len(dataset_ids), *dataset_ids]))
    # datasets interleaved per voxel, last axis fastest
    values = np.moveaxis(data, 0, -1).ravel()
    lines += [
        " ".join(f"{v:.6e}" for v in values[i : i + 6])
        for i in range(0, values.size, 6)
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


def test_read_cube_bohr_and_datasets(tmp_path):
    from molecularnodes.entities.density.cube import BOHR_TO_ANGSTROM, read_cube

    shape = (3, 4, 5)
    origin = np.array([1.0, -2.0, 0.5])
    axes = np.array([[0.4, 0.0, 0.0], [0.1, 0.5, 0.0], [0.0, 0.2, 0.6]])
    data = np.stack([_linear_field(shape), -_linear_field(shape)])
    path = _write_cube(tmp_path / "mo.cube", data, origin, axes, dataset_ids=[7, 8])

    cube = read_cube(path)
    assert cube.comments == ("comment one", "comment two")
    assert cube.dataset_ids == [7, 8]
    assert np.allclose(cube.values, data)
    assert np.allclose(cube.origin, origin * BOHR_TO_ANGSTROM)
    assert np.allclose(cube.axes, axes * BOHR_TO_ANGSTROM)
    assert cube.atomic_numbers.tolist() == [8, 1]
    assert np.allclose(cube.positions[1], [1.5 * BOHR_TO_ANGSTROM, 0, 0])


def test_read_cube_angstrom(tmp_path):
    from molecularnodes.entities.density.cube import read_cube

    shape = (2, 3, 4)
    axes = np.diag([0.3, 0.4, 0.5])
    path = _write_cube(
        tmp_path / "a.cube", _linear_field(shape)[None], [1, 2, 3], axes, angstrom=True
    )
    cube = read_cube(path)
    assert cube.values.shape == (1, *shape)
    assert np.allclose(cube.origin, [1, 2, 3])
    assert np.allclose(cube.axes, axes)
    assert np.allclose(cube.positions[1], [1.5, 0, 0])


def test_read_cube_wrong_size(tmp_path):
    from molecularnodes.entities.density.cube import read_cube

    path = _write_cube(
        tmp_path / "bad.cube", np.zeros((1, 2, 2, 2)), [0, 0, 0], np.eye(3)
    )
    path.write_text(path.read_text() + "1.0\n")
    with pytest.raises(ValueError, match="holds 9 values"):
        read_cube(path)


def test_density_affine_transform_cube(tmp_path):
    from molecularnodes.entities.density.cube import BOHR_TO_ANGSTROM

    shape = (6, 7, 8)
    origin = np.array([2.0, -1.0, 3.0])
    axes = np.array([[1.0, 0.0, 0.0], [0.5, 1.0, 0.0], [0.2, 0.3, 1.2]])
    path = _write_cube(
        tmp_path / "sheared.cube", _linear_field(shape)[None], origin, axes
    )
    density = mn.entities.density.Grids.load(path, style=None)
    origin = origin * BOHR_TO_ANGSTROM
    axes = axes * BOHR_TO_ANGSTROM
    assert np.allclose(density.axes, axes)

    # the world position of a (fractional) voxel index holds the value at that index
    ijk = np.array([[1.5, 2.25, 3.0], [4.0, 4.0, 4.0], [0.5, 6.0, 2.0], [0, 0, 0]])
    values, _ = _sample_grid_at(density, origin + ijk @ axes)
    assert np.allclose(values, ijk @ [1.0, 10.0, 100.0], atol=1e-3)


def test_cube_pyscf_aspirin(isolated_density_file):
    # aspirin HOMO from PySCF, made by docs/dev/pyscf_cubes.py
    path = isolated_density_file(data_dir / "2244_homo.cube.gz")
    density = mn.entities.density.Grids.load(path, style="density_iso_surface")
    assert density.grid.grid.shape == (24, 24, 24)
    # an orbital has both signs
    assert density.grid.grid.min() < -0.1 < 0.1 < density.grid.grid.max()

    # the atoms of the cube are the PubChem conformer they were computed from
    mol = mn.Molecule.load(path)
    ref = mn.Molecule.load(data_dir / "2244.sdf")
    assert np.allclose(
        mol.universe.atoms.positions, ref.universe.atoms.positions, atol=1e-4
    )
    assert mol.universe.atoms.elements.tolist() == ref.universe.atoms.elements.tolist()
    assert len(mol.universe.bonds) == len(ref.universe.bonds)
