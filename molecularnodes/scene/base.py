import bpy
from ..entities import Molecule
import tempfile
from pathlib import Path
import os
import warnings

try:
    from IPython.display import Image, display
except ImportError:
    Image = None
    display = None
IS_EEVEE_NEXT = bpy.app.version[0] == 4 and bpy.app.version[1] >= 2


class Canvas:
    """
    A class to handle canvas and render settings in Blender.

    This class provides properties to get and set various render settings like resolution,
    render engine, samples, FPS, and frame ranges.
    """

    def __init__(self) -> None:
        """Initialize the Canvas object."""
        pass

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
    def render_engine(self) -> str:
        """
        The current render engine for the scene.

        Returns
        -------
        str
            The name of the current render engine.
            Possible values are:
            - 'BLENDER_EEVEE': Eevee real-time render engine
            - 'CYCLES': Cycles path tracing render engine
            - 'BLENDER_WORKBENCH': Workbench render engine
        """
        return self.scene.render.engine

    @render_engine.setter
    def render_engine(self, value: str) -> None:
        """
        Set the render engine.

        Parameters
        ----------
        value : str
            The name of the render engine to set.
            Possible values are:
            - 'BLENDER_EEVEE': Eevee real-time render engine
            - 'CYCLES': Cycles path tracing render engine
            - 'BLENDER_WORKBENCH': Workbench render engine
        """

        if value == "EEVEE":
            if IS_EEVEE_NEXT:
                value = "BLENDER_EEVEE_NEXT"
            else:
                value = "BLENDER_EEVEE"

        if "EEVEE" in value and self.is_gitub_actions():
            raise ValueError("EEVEE is not supported in GitHub Actions.")

        if value == "WORKBENCH":
            value = "BLENDER_WORKBENCH"

        setattr(self.scene.render, "engine", value)

    @property
    def samples_cycles(self) -> int:
        """
        Get the number of samples for Cycles render engine.

        Returns
        -------
        int
            The number of samples for Cycles.
        """
        return self.scene.cycles.samples

    @samples_cycles.setter
    def samples_cycles(self, value: int) -> None:
        """
        Set the number of samples for Cycles render engine.

        Parameters
        ----------
        value : int
            The number of samples to set for Cycles.
        """
        self.scene.cycles.samples = value

    @property
    def samples_eevee(self) -> int:
        """
        Get the number of samples for Eevee render engine.

        Returns
        -------
        int
            The number of samples for Eevee.
        """
        return self.scene.eevee.taa_render_samples

    @samples_eevee.setter
    def samples_eevee(self, value: int) -> None:
        """
        Set the number of samples for Eevee render engine.

        Parameters
        ----------
        value : int
            The number of samples to set for Eevee.
        """
        self.scene.eevee.taa_render_samples = value

    @property
    def cycles_device(self) -> str:
        """
        Get the Cycles device type.
        """
        return self.scene.cycles.device

    @cycles_device.setter
    def cycles_device(self, value: str) -> None:
        """
        Set the Cycles device type.
        """
        value = value.upper()
        if value not in ["CPU", "GPU"]:
            raise ValueError("Invalid device type. Must be 'CPU' or 'GPU'.")
        self.scene.cycles.device = value

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

    def is_gitub_actions(self) -> bool:
        """
        Check if the script is running in GitHub Actions.
        """
        return os.getenv("GITHUB_ACTIONS") == "true"

    def frame_object(self, obj: bpy.types.Object | Molecule) -> None:
        if isinstance(obj, Molecule):
            obj = obj.object

        prev_sel = bpy.context.selected_objects
        obj.select_set(True)
        bpy.ops.view3d.camera_to_view_selected()
        obj.select_set(False)
        for o in prev_sel:
            o.select_set(True)

    def scene_reset(self) -> None:
        bpy.ops.wm.read_homefile(app_template="Molecular Nodes")

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
