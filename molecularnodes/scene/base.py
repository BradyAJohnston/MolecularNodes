import os
import tempfile
from pathlib import Path
import bpy
from ..entities import Molecule
from .engines import EEVEE, Cycles

try:
    from IPython.display import Image, display
except ImportError:
    Image = None
    display = None
IS_EEVEE_NEXT = bpy.app.version[0] == 4 and bpy.app.version[1] >= 2
IS_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"
IS_SELF_HOSTED = os.getenv("RUNNER_ENVIRONMENT") == "self-hosted"


class Canvas:
    """
    A class to handle canvas and render settings in Blender.

    This class provides properties to get and set various render settings like resolution,
    render engine, samples, FPS, and frame ranges.
    """

    def __init__(
        self,
        template: str | None = "Molecular Nodes",
        engine: EEVEE | Cycles | str = "EEVEE",
        resolution=(1280, 720),
    ) -> None:
        """Initialize the Canvas object."""
        if template:
            self.scene_reset(template=template)
        self.engine = engine
        self.resolution = resolution

    @property
    def scene(self) -> bpy.types.Scene:
        return bpy.context.scene

    @property
    def resolution(self) -> tuple[int, int]:
        """
        Get the render resolution.

        Returns
        -------
        tuple[int, int]
            A tuple containing the x and y resolution values.
        """
        x: int = self.scene.render.resolution_x
        y: int = self.scene.render.resolution_y
        return (x, y)

    @resolution.setter
    def resolution(self, value: tuple[int, int]) -> None:
        """
        Set the render resolution.

        Parameters
        ----------
        value : tuple[int, int]
            A tuple containing the x and y resolution values.
        """
        self.scene.render.resolution_x = value[0]
        self.scene.render.resolution_y = value[1]

    @property
    def engine(self) -> Cycles | EEVEE:
        return self._engine

    @engine.setter
    def engine(self, value: Cycles | EEVEE | str) -> None:
        if isinstance(value, str) and value.upper() == "CYCLES":
            self._engine = Cycles()
        elif isinstance(value, str) and value.upper() in [
            "EEVEE",
            "BLENDER_EEVEE_NEXT",
        ]:
            self._engine = EEVEE()

        else:
            self._engine = value

    @property
    def fps(self) -> float:
        """
        Get the frames per second setting.

        Returns
        -------
        float
            The current FPS value.
        """
        return self.scene.render.fps

    @fps.setter
    def fps(self, value: float) -> None:
        """
        Set the frames per second.

        Parameters
        ----------
        value : float
            The FPS value to set.
        """
        self.scene.render.fps = value

    @property
    def frame_start(self) -> int:
        """
        Get the start frame of the animation.

        Returns
        -------
        int
            The start frame number.
        """
        return self.scene.frame_start

    @frame_start.setter
    def frame_start(self, value: int) -> None:
        """
        Set the start frame of the animation.

        Parameters
        ----------
        value : int
            The start frame number to set.
        """
        self.scene.frame_start = value

    @property
    def frame_end(self) -> int:
        """
        Get the end frame of the animation.

        Returns
        -------
        int
            The end frame number.
        """
        return self.scene.frame_end

    @frame_end.setter
    def frame_end(self, value: int) -> None:
        """
        Set the end frame of the animation.

        Parameters
        ----------
        value : int
            The end frame number to set.
        """
        self.scene.frame_end = value

    def frame_object(self, obj: bpy.types.Object | Molecule) -> None:
        if isinstance(obj, Molecule):
            obj = obj.object

        prev_sel = bpy.context.selected_objects
        obj.select_set(True)
        bpy.ops.view3d.camera_to_view_selected()
        obj.select_set(False)
        for o in prev_sel:
            o.select_set(True)

    def scene_reset(
        self,
        template: str | None = "Molecular Nodes",
        engine: Cycles | EEVEE | str = "EEVEE",
    ) -> None:
        if template:
            bpy.ops.wm.read_homefile(app_template=template)
        if engine:
            self.engine = engine

    def snapshot(self, path: str | Path | None = None) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir) / "snapshot.png"

            bpy.context.scene.render.filepath = str(tmp_path)
            bpy.ops.render.render(write_still=True, animation=False)

            if display and Image:
                display(Image(tmp_path))
            else:
                # Alternative handling like saving to disk
                if path:
                    import shutil

                    shutil.copy(tmp_path, path)
