from pathlib import Path
import pytest
import molecularnodes as mn
from .constants import data_dir


def _require_image_snapshot(snapshot):
    """Return syrupy's PNG image snapshot extension or skip if unavailable."""
    try:
        from syrupy.extensions.image import PNGImageSnapshotExtension  # type: ignore
    except Exception:
        pytest.skip(
            "syrupy image extension not installed; skipping image snapshot tests"
        )
    return snapshot.use_extension(PNGImageSnapshotExtension)


def _new_canvas(resolution=(720, 480)):
    """Create a Canvas configured for fast GPU rendering; skip if no GPU."""
    canvas = mn.Canvas(engine="CYCLES", resolution=resolution)
    canvas.engine.samples = 8
    if getattr(canvas.engine, "device", "CPU") != "GPU":
        pytest.skip("GPU not available; skipping image snapshot tests")
    return canvas


MOLECULE_STYLES = [
    # Core styles
    "ball_and_stick",
    "cartoon",
    "ribbon",
    "spheres",
    "sticks",
    "surface",
    # Presets
    "preset_1",
    "preset_2",
    "preset_3",
    "preset_4",
]


@pytest.mark.parametrize("style", MOLECULE_STYLES)
def test_render_molecule_style_image(style: str, snapshot, tmp_path: Path):
    image_snapshot = _require_image_snapshot(snapshot)
    canvas = _new_canvas()

    mol = mn.Molecule.fetch("4ozs", cache=data_dir).add_style(style)
    if hasattr(mol.styles[0], "geometry"):
        mol.styles[0].geometry = "Instance"

    canvas.frame_object(mol, viewpoint="front")
    out = tmp_path / f"molecule_{style}.png"
    canvas.snapshot(path=out)

    with open(out, "rb") as f:
        assert image_snapshot == f.read()


@pytest.mark.parametrize("style", ["density_surface", "density_wire"])
def test_render_density_style_image(style: str, snapshot, tmp_path: Path):
    image_snapshot = _require_image_snapshot(snapshot)
    canvas = _new_canvas()

    density_file = data_dir / "emd_24805.map.gz"
    d1 = mn.entities.density.load(density_file, style=style)

    canvas.frame_object(d1, viewpoint="front")
    out = tmp_path / f"density_{style}.png"
    canvas.snapshot(path=out)

    with open(out, "rb") as f:
        assert image_snapshot == f.read()
