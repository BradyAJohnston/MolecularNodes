import os
import shutil
import tempfile
from contextlib import ExitStack
from enum import Enum
from pathlib import Path
from typing import Literal, Tuple
import bpy
from tqdm.auto import tqdm
from .. import assets
from ..blender import utils as blender_utils
from ..entities.base import MolecularEntity
from ..session import get_session
from ..ui import addon
from ..utils import suppress_stdout, temp_override_properties
from .camera import Camera, Viewpoints
from .engines import EEVEE, Cycles

try:
    from IPython.display import Image, Video, display
except ImportError:
    Image = None
    Video = None
    display = None

IS_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"
IS_SELF_HOSTED = os.getenv("RUNNER_ENVIRONMENT") == "self-hosted"


class ViewTransform(Enum):
    STANDARD = "Standard"
    KHRONOS = "Khronos PBR Neutral"
    AGX = "AgX"
    FILMIC = "Filmic"
    FILMIC_LOG = "Filmic Log"
    FALSE_COLOR = "False Color"
    RAW = "Raw"


_view_transform = Literal[
    ViewTransform.STANDARD.value,
    ViewTransform.KHRONOS.value,
    ViewTransform.AGX.value,
    ViewTransform.FILMIC.value,
    ViewTransform.FILMIC_LOG.value,
    ViewTransform.FALSE_COLOR.value,
    ViewTransform.RAW.value,
]


class Canvas:
    """
    A class to handle canvas and render settings in Blender.

    This class provides properties to get and set various render settings like resolution,
    render engine, samples, FPS, and frame ranges.
    """

    def __init__(
        self,
        engine: EEVEE | Cycles | str = "EEVEE",
        resolution=(1280, 720),
        transparent: bool = False,
        template: str | None = "Molecular Nodes",
    ) -> None:
        """Initialize the Canvas object."""
        addon.register()
        assets.install()
        if template:
            self.scene_reset(template=template)
        self.engine = engine
        self.resolution = resolution
        self.camera = Camera()
        self.transparent = transparent

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
        if isinstance(value, Cycles) or isinstance(value, EEVEE):
            self._engine = value
        elif isinstance(value, str):
            if value.upper() == "CYCLES":
                self._engine = Cycles()
            elif value.upper() in [
                "EEVEE",
                "BLENDER_EEVEE_NEXT",
                "BLENDER_EEVEE",
            ]:
                self._engine = EEVEE()
            else:
                raise ValueError("String does not match either 'EEVEE' or 'CYCLES'")
        else:
            raise ValueError(
                "Must be either a string selecting the render engine or a dataclass mn.scene.Cycles()"
            )

        self._engine._enable_engine()

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

    @property
    def transparent(self) -> bool:
        """
        Get the transparency setting for rendering.

        Returns
        -------
        bool
            True if transparency is enabled, False otherwise.
        """
        return self.scene.render.film_transparent

    @transparent.setter
    def transparent(self, value: bool) -> None:
        """
        Set the transparency setting for rendering.

        Parameters
        ----------
        value : bool
            True to enable transparency, False to disable.
        """
        self.scene.render.film_transparent = value

    @property
    def background(self) -> Tuple[float, float, float, float]:
        return (
            self.scene.world.node_tree.nodes["MN_world_shader"].inputs[3].default_value
        )

    @background.setter
    def background(self, value: Tuple[float, float, float, float]) -> None:
        self.scene.world.node_tree.nodes["MN_world_shader"].inputs[
            3
        ].default_value = value

    # @property
    # def hdri_strength

    @property
    def view_transform(self) -> _view_transform:
        """
        Get the current view transform setting.

        Returns
        -------
        str
            The current view transform value.
        """
        return self.scene.view_settings.view_transform

    @view_transform.setter
    def view_transform(self, value: _view_transform | ViewTransform | str) -> None:
        """
        Set the view transform setting.

        Parameters
        ----------
        value : str | ViewTransform
            The view transform value to set. Accepts enum, full name, or lowercase/shortened name.
        """
        # Normalize input
        if isinstance(value, ViewTransform):
            vt_value = value.value
        elif isinstance(value, str):
            value_lower = value.strip().lower()
            # Try to match full name (case-insensitive)
            for vt in ViewTransform:
                if value_lower == vt.value.lower():
                    vt_value = vt.value
                    break
            else:
                # Try to match shortened name (e.g., "standard", "agx", "filmic", etc.)
                for vt in ViewTransform:
                    if value_lower == vt.name.lower():
                        vt_value = vt.value
                        break
                else:
                    raise ValueError(
                        f"Invalid view transform '{value}'. Must be one of {[vt.value for vt in ViewTransform]} or a valid enum name."
                    )
        else:
            raise TypeError("view_transform must be a str or ViewTransform enum.")

        self.scene.view_settings.view_transform = vt_value

    def frame_object(
        self, obj: bpy.types.Object | MolecularEntity, viewpoint: Viewpoints = None
    ) -> None:
        """
        Frame an object or Molecular entity

        Parameters
        ----------
        obj : bpy.types.Object | MolecularEntity
            Blender object or Molecular entity to frame.

        viewpoint: str, optional
            Viewing direction along an axis
            One of ["default", "front", "back", "top", "bottom", "left", "right"]

        """
        if isinstance(obj, MolecularEntity):
            obj = obj.object
        # set the camera viewpoint if specified
        if viewpoint is not None:
            self.camera.set_viewpoint(viewpoint)
        # set the camera to look at the object
        blender_utils.look_at_object(obj)

    def frame_view(
        self, view: list[tuple] | MolecularEntity, viewpoint: Viewpoints = None
    ) -> None:
        """
        Frame one or more views of Molecular entities
        Multiple views can be added with + to combine into a single view

        Parameters
        ----------
        view: list[tuple]
            A bounding box (set of 8 3D vertices) of the region of interest

        viewpoint: str, optional
            Viewing direction along an axis
            One of ["default", "front", "back", "top", "bottom", "left", "right"]

        """
        if isinstance(view, MolecularEntity):
            view_tuple = view.get_view()
        else:
            view_tuple = view
        # set the camera viewpoint if specified
        if viewpoint is not None:
            self.camera.set_viewpoint(viewpoint)
        # set the camera to look at the bounding box of the view
        blender_utils.look_at_bbox(view_tuple)

    def clear(self) -> None:
        """
        Clear all entities from the scene without affecting lighting or render settings.
        """
        session = get_session()
        # Iterate over a list copy to avoid modifying the dict during iteration
        for entity in list(session.entities.values()):
            session.remove(entity.uuid)

    def scene_reset(
        self,
        template: str | None = "Molecular Nodes",
        engine: Cycles | EEVEE | str = "EEVEE",
    ) -> None:
        if template:
            bpy.ops.wm.read_homefile(app_template=template)
        if engine:
            self.engine = engine

    def snapshot(
        self,
        path: str | Path | None = None,
        frame: int | None = None,
        file_format: str = "PNG",
    ) -> None:
        """
        Render an image of the current scene.

        Parameters
        ----------
        path : str | Path | None, optional
            File path to write the rendered image to.

        frame : int, optional
            Frame number of scene to render. When not specified,
            current scene's current_frame is used

        file_format : str, optional
            File format of the rendered image.

        """
        scene = self.scene
        render_settings = scene.render
        image_settings = render_settings.image_settings
        render_frame = scene.frame_current if frame is None else frame
        with ExitStack() as stack:
            # set the use_file_extension to auto generate file extension
            # set the file_format to the specified one
            temp_override_properties(
                stack,
                [
                    (render_settings, "use_file_extension", True),
                    (image_settings, "file_format", file_format),
                    (scene, "frame_current", render_frame),
                ],
            )
            # create temporary file with the file_format extension
            tmp_file = stack.enter_context(
                tempfile.NamedTemporaryFile(suffix=render_settings.file_extension)
            )
            # set the filepath to the temporary file
            temp_override_properties(
                stack,
                [
                    (render_settings, "filepath", tmp_file.name),
                ],
            )
            # set the frame number - only frame_set will update animation data
            # scene.current_frame will only update the timeline
            scene.frame_set(render_frame)
            # suppress stdout output generated by render process
            with suppress_stdout():
                # render the image
                bpy.ops.render.render(write_still=True, animation=False)
            if path:
                # save to file if path specified
                shutil.copy(tmp_file.name, path)
            elif display and Image:
                # only display in notebook if path not specified
                # and in notebook context
                display(Image(tmp_file.name))

    def animation(
        self,
        path: str | Path | None = None,
        frame_start: int | None = None,
        frame_end: int | None = None,
        render_scale: int = 100,
    ) -> None:
        """
        Render an animation of the current scene.

        Parameters
        ----------
        path : str | Path | None, optional
            File path to write the rendered animation to.

        frame_start : int, optional
            Start frame of the animation. When not specified, current scene's
            start frame is used

        frame_end : int, optional
            End frame of the animation. When not specified, current scene's
            end frame is used

        render_scale : int, optional
            Scale of the rendered animation frames with respect to the resolution.

        """
        # determine frame range
        start = self.frame_start if frame_start is None else frame_start
        end = self.frame_end if frame_end is None else frame_end
        if end < start:
            raise ValueError(f"End frame {end} cannot be less than start frame {start}")
        n = len(str(end))
        frame_range = range(start, end + 1)

        scene = self.scene
        render_settings = scene.render
        image_settings = render_settings.image_settings
        # create a temporary directory
        with tempfile.TemporaryDirectory() as tmp_dir:
            # render individual frames
            with ExitStack() as stack:
                temp_override_properties(
                    stack,
                    [
                        (render_settings, "use_lock_interface", True),
                        (render_settings, "resolution_percentage", render_scale),
                        (render_settings, "filepath", ""),
                        (image_settings, "file_format", "PNG"),
                        (scene, "frame_current", start),
                    ],
                )
                it = tqdm(frame_range, desc="Rendering frames")
                # render individual frames
                for frame in it:
                    render_image = os.path.join(tmp_dir, str(frame).zfill(n)) + ".png"
                    # set the output file
                    scene.render.filepath = render_image
                    # set the frame - only frame_set updates animation data
                    scene.frame_set(frame)
                    with suppress_stdout():
                        # render each frame image
                        bpy.ops.render.render(write_still=True)
                it.close()

            # add to video sequence editor
            sequence_editor = scene.sequence_editor_create()
            strips = sequence_editor.strips  # .sequences is deprecated
            image_strip = None
            for frame in frame_range:
                render_image = os.path.join(tmp_dir, str(frame).zfill(n)) + ".png"
                if image_strip is None:
                    image_strip = strips.new_image("name", render_image, 1, start)
                else:
                    image_strip.elements.append(os.path.basename(render_image))

            # render animation
            video_file = os.path.join(tmp_dir, "animation.mp4")
            with ExitStack() as stack:
                temp_override_properties(
                    stack,
                    [
                        (render_settings, "use_lock_interface", True),
                        (render_settings, "resolution_percentage", render_scale),
                        (render_settings, "filepath", ""),
                        (image_settings, "file_format", "FFMPEG"),
                        (render_settings.ffmpeg, "format", "MPEG4"),
                        (scene, "frame_current", start),
                        (scene, "frame_start", start),
                        (scene, "frame_end", end),
                    ],
                )
                scene.render.filepath = video_file
                it = tqdm(range(0, 1), desc="Generating video")
                for i in it:
                    with suppress_stdout():
                        # render the animation of sequence editor images
                        bpy.ops.render.render(animation=True)
                it.close()

            # clear the video sequence editor
            strips.remove(image_strip)

            if path:
                # save to file if path specified
                shutil.copy(video_file, path)
            elif display and Video:
                # only display in notebook if path not specified
                # and in notebook context
                display(Video(video_file, embed=True))
