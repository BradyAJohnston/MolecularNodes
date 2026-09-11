"""
Golden-image tests for final renders.

Renders are compared against committed reference PNGs with the same tolerance
model Blender's own render tests use (see `compare_images` in utils.py): a
pixel fails when a channel differs by more than ~4/255, and the image fails
when more than 1% of pixels fail. This absorbs the small cross-platform
differences of CPU rendering while still catching real regressions in styles,
materials, lighting, color management and compositing.

Update the reference images with `pytest --snapshot-update` after an
intentional visual change, and eyeball the new goldens before committing. On
failure the received and difference images are written to
``tests/image_failures/``.
"""

from itertools import product
import pytest
import molecularnodes as mn
import molecularnodes.nodes.geometry as mg
from .constants import data_dir
from .utils import ImageSnapshotExtension


@pytest.fixture
def image_snapshot(snapshot):
    return snapshot.use_extension(ImageSnapshotExtension)


@pytest.fixture
def golden_canvas():
    # deterministic render settings: Cycles on the CPU with a fixed seed and
    # sample count, and no denoising (OpenImageDenoise output is not stable
    # across versions and platforms), compositor on the CPU for GPU-less CI
    # runners. Same-seed CPU renders are pixel-identical per platform, so the
    # comparison tolerance is entirely headroom for cross-platform drift; the
    # generous sample count keeps residual noise well below the threshold.
    canvas = mn.Canvas(resolution=(128, 128))
    canvas.engine = mn.scene.Cycles(samples=256, device="CPU", denoise=False)
    canvas.compositor.device = "CPU"
    return canvas


def _render(canvas, tmp_path) -> bytes:
    file = tmp_path / "render.png"
    canvas.snapshot(file)
    return file.read_bytes()


def _fetch_molecule():
    return mn.Molecule.fetch("4ozs", cache=data_dir, format="bcif")


def test_render_spheres(golden_canvas, tmp_path, image_snapshot):
    mol = _fetch_molecule()
    mol.add_style("spheres")
    golden_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(golden_canvas, tmp_path)


def test_render_cartoon(golden_canvas, tmp_path, image_snapshot):
    mol = _fetch_molecule()
    mol.add_style("cartoon")
    golden_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(golden_canvas, tmp_path)


def test_render_annotations(golden_canvas, tmp_path, image_snapshot):
    # exercises the annotation overlay through the compositor
    mol = _fetch_molecule()
    mol.add_style("cartoon")
    mol.annotations.add_atom_info(selection="name CA and resid 20", show_resid=True)
    golden_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(golden_canvas, tmp_path)


# spheres and cartoon are covered by the dedicated tests above
@pytest.mark.parametrize("style", ["ball_and_stick", "ribbon", "sticks", "surface"])
def test_render_style(style, golden_canvas, tmp_path, image_snapshot):
    mol = _fetch_molecule()
    mol.add_style(style)
    golden_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(golden_canvas, tmp_path)


@pytest.mark.parametrize(
    "material",
    ["Default", "AmbientOcclusion", "Flat", "Squishy", "TransparentOutline"],
)
def test_render_material(material, golden_canvas, tmp_path, image_snapshot):
    mol = _fetch_molecule()
    mol.add_style("surface", material=getattr(mn.material, material)())
    golden_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(golden_canvas, tmp_path)


def test_render_selection_string(golden_canvas, tmp_path, image_snapshot):
    # a selection phrase masks the second style to part of the molecule
    mol = _fetch_molecule()
    mol.add_style("cartoon")
    mol.add_style("spheres", selection="resid 1:40")
    golden_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(golden_canvas, tmp_path)


def test_render_selection_atomgroup(golden_canvas, tmp_path, image_snapshot):
    # an MDAnalysis AtomGroup can be used as a selection directly
    mol = _fetch_molecule()
    mol.add_style("cartoon")
    mol.add_style("sticks", selection=mol.universe.select_atoms("resid 100:150"))
    golden_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(golden_canvas, tmp_path)


@pytest.fixture
def assembly_image_snapshot(snapshot):
    # full-frame assemblies (the 8OZK capsid especially) are almost entirely
    # sub-pixel ribbon geometry at 128x128, so cross-platform floating-point
    # drift flips pixel coverage across the whole silhouette: up to ~10% of
    # pixels differ between platforms while the per-channel error stays small
    # (diff images are uniform speckle, no structural difference). Allow a
    # higher failing-pixel percentage here; a real assembly regression (wrong
    # scale, missing instances) changes large regions and still fails.
    return snapshot.use_extension(
        ImageSnapshotExtension.with_tolerance(fail_percent=15)
    )


@pytest.mark.parametrize(
    "code,node,assembly",
    list(
        product(
            ["4ozs", "1cd3", "8OZK"], ["code", "method"], ["assembly", "no_assembly"]
        )
    ),
)
def test_render_assembly(
    golden_canvas, tmp_path, assembly_image_snapshot, code, node, assembly
):
    mol = mn.Molecule.fetch(code)
    mat = mn.material.Flat()
    if node == "code":
        with mol.tree.reset() as (atoms, join):
            (
                atoms
                >> mg.StyleRibbon(material=mat.material)
                >> (
                    mg.AssemblyInstance(data_object=mol.create_data_object())
                    if assembly == "assembly"
                    else None
                )
                >> join
            )
    else:
        mol.add_style("ribbon", material=mat, assembly=(assembly == "assembly"))
    golden_canvas.look_at(mol, viewpoint="front")
    assert assembly_image_snapshot == _render(golden_canvas, tmp_path)


@pytest.mark.parametrize(
    "files,backbone",
    list(
        product(
            (
                ["minicircle.top", "minicircle.dat"],
                ["linear.top", "linear_traj.dat"],
                ["origami_old.top", "origami_old.dat"],
            ),
            ("Arrows", "Curve"),
        )
    ),
)
def test_render_oxdna_simple_circle(
    golden_canvas, tmp_path, assembly_image_snapshot, files, backbone
):
    ens = mn.entities.OXDNA.load(
        topology=data_dir / f"oxdna/{files[0]}",
        coordinates=data_dir / f"oxdna/{files[1]}",
    )

    with ens.tree.reset() as (atoms, join):
        (
            atoms
            >> mg.CentreOnSelection()
            >> mg.OxDNAStyleRibbon(
                backbone_shape=backbone, material=mn.material.Default().material
            )
            >> join
        )

    golden_canvas.look_at(ens, viewpoint="top")
    assert assembly_image_snapshot == _render(golden_canvas, tmp_path)


@pytest.mark.parametrize(
    "pdb_id,style",
    list(product(("1KAL",), ("ribbon", "ball_and_stick", "cartoon"))),
)
def test_render_cyclic_peptide(
    golden_canvas, tmp_path, assembly_image_snapshot, pdb_id, style
):
    mol = mn.Molecule.fetch(pdb_id).add_style(style)
    golden_canvas.look_at(mol, viewpoint="top")
    assert assembly_image_snapshot == _render(golden_canvas, tmp_path)
