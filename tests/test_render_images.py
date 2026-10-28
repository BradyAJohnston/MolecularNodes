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
import bpy
import pytest
from nodebpy.nodes.geometry import RandomValue, StoreNamedAttribute
import molecularnodes as mn
import molecularnodes.nodes.compositor as mc
import molecularnodes.nodes.geometry as mg
from molecularnodes.material import (
    AmbientOcclusion,
    Default,
    Flat,
    Squishy,
    Transparent,
)
from .constants import data_dir
from .test_nodes import SYMMETRY_EXAMPLES
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


def _fetch_molecule(fetch):
    return fetch("4ozs", format="bcif")


def test_render_spheres(fetch, golden_canvas, tmp_path, image_snapshot):
    mol = _fetch_molecule(fetch)
    mol.add_style("spheres")
    golden_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(golden_canvas, tmp_path)


def test_render_cartoon(fetch, golden_canvas, tmp_path, image_snapshot):
    mol = _fetch_molecule(fetch)
    mol.add_style("cartoon")
    golden_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(golden_canvas, tmp_path)


def test_render_annotations(fetch, golden_canvas, tmp_path, image_snapshot):
    # exercises the annotation overlay through the compositor
    mol = _fetch_molecule(fetch)
    mol.add_style("cartoon")
    mol.annotations.add_atom_info(selection="name CA and resid 20", show_resid=True)
    golden_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(golden_canvas, tmp_path)


# spheres and cartoon are covered by the dedicated tests above
@pytest.mark.parametrize("style", ["ball_and_stick", "ribbon", "sticks", "surface"])
def test_render_style(fetch, style, golden_canvas, tmp_path, image_snapshot):
    mol = _fetch_molecule(fetch)
    mol.add_style(style)
    golden_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(golden_canvas, tmp_path)


@pytest.mark.parametrize(
    "material",
    ["Default", "AmbientOcclusion", "Flat", "Squishy", "Transparent"],
)
def test_render_material(fetch, material, golden_canvas, tmp_path, image_snapshot):
    mol = _fetch_molecule(fetch)
    mol.add_style("surface", material=getattr(mn.material, material)())
    golden_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(golden_canvas, tmp_path)


def test_render_material_flat_outline(fetch, golden_canvas, tmp_path, image_snapshot):
    # Flat has no outline by default; outline=True adds the dark edge lines
    mol = _fetch_molecule(fetch)
    mol.add_style("surface", material=mn.material.Flat(outline=True))
    golden_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(golden_canvas, tmp_path)


def test_render_selection_string(fetch, golden_canvas, tmp_path, image_snapshot):
    # a selection phrase masks the second style to part of the molecule
    mol = _fetch_molecule(fetch)
    mol.add_style("cartoon")
    mol.add_style("spheres", selection="resid 1:40")
    golden_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(golden_canvas, tmp_path)


def test_render_selection_atomgroup(fetch, golden_canvas, tmp_path, image_snapshot):
    # an MDAnalysis AtomGroup can be used as a selection directly
    mol = _fetch_molecule(fetch)
    mol.add_style("cartoon")
    mol.add_style("sticks", selection=mol.universe.select_atoms("resid 100:150"))
    golden_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(golden_canvas, tmp_path)


def _outline(canvas, **inputs):
    from nodebpy import compositor as c

    canvas.passes = ["combined", "z", "normal"]
    with canvas.compositor.reset() as (image, output):
        render = c.RenderLayers()
        mask = mc.DepthOutline(
            depth=render.o.depth, normals_pass=render.o.normal, **inputs
        )
        mc.Outline(image=image, mask=mask) >> output


def _chain_outline(canvas):
    from nodebpy import compositor as c

    canvas.passes = ["combined"]
    canvas.add_aov("chain_id")
    for mat in bpy.data.materials:
        if mat.name.startswith("Flat"):
            mn.material.add_aov(mat, "chain_id", offset=1.0)
    with canvas.compositor.reset() as (image, output):
        render = c.RenderLayers()
        mask = mc.IDOutline(id=render.o.chain_id)
        mc.Outline(image=image, mask=mask) >> output


def _fog(canvas, image, depth, output):
    near, far = mn.scene.compositor._depth_range(canvas.scene)
    (
        mc.CompositeDepthFog(
            image=image,
            depth=depth,
            near=near,
            far=far,
            back_fog=0.2,
            fog_color=(0.2, 0.3, 0.8, 1.0),
        )
        >> output
    )


def _fog_only(canvas):
    from nodebpy import compositor as c

    canvas.passes = ["combined", "z"]
    with canvas.compositor.reset() as (image, output):
        render = c.RenderLayers()
        _fog(canvas, image, render.o.depth, output)


def _illustrate_fog(canvas):
    from nodebpy import compositor as c

    canvas.passes = ["combined", "z"]
    with canvas.compositor.reset() as (image, output):
        render = c.RenderLayers()
        node = mc.Illustrate(
            image=image,
            alpha=render.o.alpha,
            depth=render.o.depth,
            shadow=mc.ConeShadow(
                depth=render.o.depth,
                pixel_size=mn.scene.compositor._pixel_size_angstrom(canvas.scene)
                * mn.scene.compositor.WORLD_SCALE,
            ),
        )
        _fog(canvas, node, render.o.depth, output)


# name -> (material, compositor setup applied after the camera is framed)
COMPOSITOR_CASES = {
    "plain": ("Default", lambda canvas: None),
    # the two Depth Outline methods, then the normal edges and dilation
    "outline": ("Flat", _outline),
    "outline_goodsell": ("Flat", lambda canvas: _outline(canvas, method="Goodsell")),
    "outline_normals": (
        "Default",
        lambda canvas: _outline(canvas, normals=True, size=2),
    ),
    "chain_outline": ("Flat", _chain_outline),
    "fog": ("Default", _fog_only),
    "illustrate_cone": ("Flat", lambda canvas: canvas.compositor.illustrate()),
    "illustrate_goodsell": (
        "Flat",
        lambda canvas: canvas.compositor.illustrate(outline_method="Goodsell"),
    ),
    "illustrate_ao": (
        "Default",
        lambda canvas: canvas.compositor.illustrate(shading="ao", outline_normals=True),
    ),
    "illustrate_fog": ("Flat", _illustrate_fog),
}


# cases rendered on another structure than 4ozs: the bases of 1BNA's cartoon are
# instances, so its chain outline only separates the two strands when the AOV is
# read from the instancer
CASE_STRUCTURES = {"chain_outline": "1BNA"}


@pytest.mark.parametrize("case", list(COMPOSITOR_CASES))
def test_render_compositor(case, fetch, golden_canvas, tmp_path, image_snapshot):
    # the compositor node assets: depth and chain outlines, depth
    # fog and the Illustrate composite with each of its shading sources
    material, setup = COMPOSITOR_CASES[case]
    golden_canvas.transparent = True
    mol = fetch(CASE_STRUCTURES.get(case, "4ozs"), format="bcif")
    mol.add_style("cartoon", material=getattr(mn.material, material)())
    golden_canvas.look_at(mol, viewpoint="front")
    setup(golden_canvas)
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
    fetch, golden_canvas, tmp_path, assembly_image_snapshot, code, node, assembly
):
    mol = fetch(code)
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
            ("Sticks", "Ribbon"),
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
            >> mg.OxDNAStyleClassic(
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
    fetch, golden_canvas, tmp_path, assembly_image_snapshot, pdb_id, style
):
    mol = fetch(pdb_id).add_style(style)
    golden_canvas.look_at(mol, viewpoint="top")
    assert assembly_image_snapshot == _render(golden_canvas, tmp_path)


@pytest.mark.parametrize("code", list(SYMMETRY_EXAMPLES))
def test_render_symmetry(fetch, golden_canvas, tmp_path, assembly_image_snapshot, code):
    # generated symmetry reproducing the entry's deposited assembly; the rings
    # are viewed down their axis and the helix side on
    mol = fetch(code)
    with mol.tree.reset() as (atoms, join):
        (
            atoms
            >> mg.StyleRibbon(material=mn.material.Flat().material)
            >> SYMMETRY_EXAMPLES[code]()
            >> StoreNamedAttribute.instance.color(
                name="Color",
                value=mg.ColorMatplotlib(
                    RandomValue(), category="Uniform", uniform="magma"
                ),
            )
            >> StoreNamedAttribute.instance.boolean(name="is_instanced", value=True)
            >> join
        )
    golden_canvas.look_at(mol, viewpoint="front" if code == "4UDV" else "top")
    assert assembly_image_snapshot == _render(golden_canvas, tmp_path)


@pytest.mark.parametrize(
    "material",
    [mn.material.Default, mn.material.AmbientOcclusion],
)
def test_render_faded_geometry(
    fetch, golden_canvas, tmp_path, assembly_image_snapshot, material
):
    mol = fetch("1een")
    with mol.tree.reset() as (atoms, join):
        (atoms >> mg.StyleCartoon(material=mn.material.Flat().material) >> join)

        (
            atoms
            >> mg.StyleSurface(material=material().material)
            >> mg.FadeGeometry(fade=0.3)
            >> join
        )

    golden_canvas.look_at(mol)
    assert assembly_image_snapshot == _render(golden_canvas, tmp_path)


@pytest.mark.parametrize(
    "viewpoint",
    ["front", "pca", "unobstructed", "visible_area", "raycast", "entropy"],
)
def test_render_viewpoint(fetch, golden_canvas, tmp_path, image_snapshot, viewpoint):
    # retinoic acid sits inside a beta barrel: the front view sees it through
    # a strand, the occlusion-aware viewpoints should look in through the
    # opening; entropy picks a direction for the whole protein
    mol = fetch("1cbs")
    mol.add_style("cartoon", selection="protein", color=(0.6, 0.6, 0.65, 1.0))
    mol.add_style("ball_and_stick", selection="resname REA")
    ligand = mol.get_view("resname REA")
    golden_canvas.look_at(
        mol if viewpoint == "entropy" else ligand, viewpoint=viewpoint
    )
    # keep the chosen direction, frame the ligand and its pocket
    golden_canvas.look_at(
        mol.get_view("resname REA or (protein and around 6 resname REA)")
    )
    assert image_snapshot == _render(golden_canvas, tmp_path)


def _style_surface(m):
    return lambda: mn.nodes.geometry.StyleSurface(quality=4, scale=1.0, material=m)


def _material_flat():
    return mn.material.Flat(outline=True)


@pytest.mark.parametrize(
    "name,style,material",
    list(
        product(
            ["caffeine"],
            [
                "spheres",
                "sticks",
                _style_surface,
                "ball_and_stick",
            ],
            [
                mn.material.Default,
                _material_flat,
                mn.material.AmbientOcclusion,
            ],
        )
    ),
)
def test_render_small_molecule(
    fetch, golden_canvas, tmp_path, image_snapshot, name, style, material
):
    mol = fetch(name, format="sdf", database="pubchem")
    if isinstance(style, str):
        mol.add_style(style, material=material().material)
    else:
        mol.add_style(style(material().material))
    golden_canvas.look_at(mol, viewpoint="top")
    assert image_snapshot == _render(golden_canvas, tmp_path)


@pytest.fixture
def eevee_canvas():
    # EEVEE needs a GPU, so tests using this are marked `eevee` and skipped on
    # machines without one (most CI runners). A fixed sample count keeps the
    # renders deterministic on a given GPU.
    canvas = mn.Canvas(resolution=(128, 128))
    canvas.engine = mn.scene.EEVEE(samples=64)
    canvas.compositor.device = "CPU"
    return canvas


@pytest.mark.eevee
@pytest.mark.parametrize(
    "style", ["spheres", "cartoon", "ball_and_stick", "ribbon", "sticks", "surface"]
)
def test_render_eevee_style(fetch, style, eevee_canvas, tmp_path, image_snapshot):
    # spheres are instanced meshes under EEVEE rather than a point cloud
    mol = _fetch_molecule(fetch)
    mol.add_style(style)
    eevee_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(eevee_canvas, tmp_path)


@pytest.mark.eevee
@pytest.mark.parametrize(
    "material",
    [Default, AmbientOcclusion, Flat, Squishy, Transparent],
)
def test_render_eevee_material(fetch, material, eevee_canvas, tmp_path, image_snapshot):
    mol = _fetch_molecule(fetch)
    mol.add_style("surface", material=material())
    eevee_canvas.look_at(mol, viewpoint="front")
    assert image_snapshot == _render(eevee_canvas, tmp_path)


@pytest.mark.eevee
def test_render_eevee_illustrate_shadow(fetch, eevee_canvas, tmp_path, image_snapshot):
    # the shadow pass only exists under EEVEE, Cycles can't cover this
    eevee_canvas.transparent = True
    mol = _fetch_molecule(fetch)
    mol.add_style("cartoon", material=mn.material.Flat())
    eevee_canvas.look_at(mol, viewpoint="front")
    eevee_canvas.compositor.illustrate(shading="shadow")
    assert image_snapshot == _render(eevee_canvas, tmp_path)
