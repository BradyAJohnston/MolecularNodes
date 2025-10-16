"""
Image snapshot tests for visual regression testing.

Use assert_image_snapshot() to compare rendered images with stored snapshots.

Example:
    @pytest.mark.parametrize("style", ["ball_and_stick", "cartoon"])
    def test_molecule_render(style, image_snapshot, tmp_path, request):
        canvas = _new_canvas()
        mol = mn.Molecule.fetch("4ozs").add_style(style)

        # Render
        canvas.frame_object(mol, viewpoint="front")
        out = tmp_path / f"{style}.png"
        canvas.snapshot(path=out)

        # Compare with snapshot - handles everything automatically
        assert_image_snapshot(out, image_snapshot, request)
"""
from pathlib import Path
import pytest
from PIL import Image
import molecularnodes as mn
from .constants import data_dir

# Mark all tests in this file as image_snapshot tests
pytestmark = pytest.mark.image_snapshot


def _new_canvas(resolution=(720, 480)):
    """Create a Canvas configured for fast GPU rendering; skip if no GPU."""
    canvas = mn.Canvas(engine="CYCLES", resolution=resolution)
    canvas.engine.samples = 8
    if getattr(canvas.engine, "device", "CPU") != "GPU":
        pytest.skip("GPU not available; skipping image snapshot tests")
    return canvas


def assert_image_snapshot(image_or_path, image_snapshot, request, threshold=0.01):
    """
    Compare a rendered image with a stored snapshot.

    Handles:
    - Automatic snapshot path generation from test name (including parameters)
    - Loading image from file path or using existing PIL Image
    - Calling image_snapshot with correct arguments

    Args:
        image_or_path: Either a file path (str/Path) or PIL Image object
        image_snapshot: The pytest image_snapshot fixture
        request: The pytest request fixture
        threshold: Pixelmatch threshold for comparison (default 0.01)

    Example with file path:
        canvas.snapshot(path=out)
        assert_image_snapshot(out, image_snapshot, request)

    Example with PIL Image:
        img = Image.open(path)
        assert_image_snapshot(img, image_snapshot, request, threshold=0.02)
    """
    # Load image if a path was provided
    if isinstance(image_or_path, (str, Path)):
        image = Image.open(image_or_path)
    else:
        image = image_or_path

    # Generate snapshot path from test name (includes parameters)
    test_name = request.node.name
    snapshot_path = f"tests/__image_snapshots__/{test_name}.png"

    # Compare with snapshot
    image_snapshot(image, snapshot_path, threshold=threshold)


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
def test_render_molecule_style_image(style: str, image_snapshot, tmp_path: Path, request):
    """Test molecule rendering with different styles."""
    canvas = _new_canvas()

    mol = mn.Molecule.fetch("4ozs", cache=data_dir).add_style(style)
    if hasattr(mol.styles[0], "geometry"):
        mol.styles[0].geometry = "Instance"

    canvas.frame_object(mol, viewpoint="front")
    out = tmp_path / f"molecule_{style}.png"
    canvas.snapshot(path=out)

    assert_image_snapshot(out, image_snapshot, request)


@pytest.mark.parametrize("style", ["density_surface", "density_wire"])
def test_render_density_style_image(style: str, image_snapshot, tmp_path: Path, request):
    """Test density rendering with different styles."""
    canvas = _new_canvas()

    density_file = data_dir / "emd_24805.map.gz"
    d1 = mn.entities.density.load(density_file, style=style)

    canvas.frame_object(d1, viewpoint="front")
    out = tmp_path / f"density_{style}.png"
    canvas.snapshot(path=out)

    assert_image_snapshot(out, image_snapshot, request)
