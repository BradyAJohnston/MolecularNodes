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
from ..assets.template import list_templates
from ..blender import utils as blender_utils
from ..entities.base import MolecularEntity
from ..scene.compositor import setup_compositor
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
    High-level render controller for Blender scenes.

    Canvas configures the active Blender scene for Molecular Nodes renders
    (engine, resolution, transparency, color management), exposes convenient
    properties for common render settings, and provides helpers to frame
    objects/views and render stills or animations.

    Parameters
    ----------
    engine : EEVEE | Cycles | str, default "EEVEE"
        Render engine to use. Accepts an instance of ``mn.scene.EEVEE`` or
        ``mn.scene.Cycles``, or a case-insensitive string: ``"EEVEE"`` or
        ``"CYCLES"``.
    resolution : tuple[int, int], default (1280, 720)
        Output resolution in pixels as ``(width, height)``.
    transparent : bool, default False
        When ``True``, renders use a transparent film (alpha background).
    template : pathlib.Path | str | None, default "Molecular Nodes"
        Scene template to initialize. If a string is provided it can be either
        the name of an installed Blender app template (e.g. ``"Molecular Nodes"``),
        or a path to a ``.blend`` file. If ``None``, the Blender default startup
        file is used.

    Attributes
    ----------
    scene : bpy.types.Scene
        The active Blender scene.
    camera : molecularnodes.scene.camera.Camera
        Convenience camera controller bound to the active scene camera.
    engine : EEVEE | Cycles
        The configured render engine object.
    resolution : tuple[int, int]
        Current render resolution in pixels.
    transparent : bool
        Whether the film background is transparent.
    fps : float
        Frames per second for animation output.
    frame_start : int
        Start frame of the scene range.
    frame_end : int
        End frame of the scene range.
    background : tuple[float, float, float, float]
        World background color as RGBA in the range [0, 1].
    view_transform : {"Standard", "Khronos PBR Neutral", "AgX", "Filmic", "Filmic Log", "False Color", "Raw"}
        Active view transform for color management.

    Examples
    --------
    Create a canvas and render a snapshot with a transparent background::

        import molecularnodes as mn
        cv = mn.Canvas(engine="CYCLES", resolution=(800, 800), transparent=True)
        cv.snapshot("frame.png")

    See Also
    --------
    molecularnodes.scene.engines.EEVEE : Render engine configuration for EEVEE.
    molecularnodes.scene.engines.Cycles : Render engine configuration for Cycles.
    molecularnodes.scene.camera.Camera : Camera controller used by Canvas.
    """

    def __init__(
        self,
        engine: EEVEE | Cycles | str = "EEVEE",
        resolution=(1280, 720),
        transparent: bool = False,
        template: Path | str | None = "Molecular Nodes",
    ) -> None:
        """
        Initialize a Canvas and prepare the Blender scene.

        Parameters
        ----------
        engine : EEVEE | Cycles | str, default "EEVEE"
            Render engine configuration or engine name.
        resolution : tuple[int, int], default (1280, 720)
            Output resolution in pixels as ``(width, height)``.
        transparent : bool, default False
            Enable a transparent film (alpha background) when ``True``.
        template : pathlib.Path | str | None, default "Molecular Nodes"
            Scene template name or path to a ``.blend`` file. Use ``None`` for
            Blender's default startup scene.
        """
        addon.register()
        assets.install()
        if template:
            self.scene_reset(template=template)
        self.engine = engine
        self.resolution = resolution
        self.camera = Camera()
        self.transparent = transparent
        setup_compositor(self.scene)

    @property
    def scene(self) -> bpy.types.Scene:
        """
        Get the active Blender scene.

        Returns
        -------
        bpy.types.Scene
            The current context scene.
        """
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
        """
        Get the configured render engine.

        Returns
        -------
        EEVEE | Cycles
            The active engine configuration object.
        """
        return self._engine

    @engine.setter
    def engine(self, value: Cycles | EEVEE | str) -> None:
        """
        Set the render engine.

        Parameters
        ----------
        value : EEVEE | Cycles | str
            Either an engine configuration instance or a case-insensitive
            string: ``"EEVEE"`` or ``"CYCLES"``.

        Raises
        ------
        ValueError
            If an unsupported string is provided.
        ValueError
            If the value is neither a valid string nor supported engine type.
        """
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
        """
        Get the world background color.

        Returns
        -------
        tuple[float, float, float, float]
            RGBA values in the range [0, 1].
        """
        return (
            self.scene.world.node_tree.nodes["MN_world_shader"].inputs[3].default_value
        )

    @background.setter
    def background(self, value: Tuple[float, float, float, float]) -> None:
        """
        Set the world background color.

        Parameters
        ----------
        value : tuple[float, float, float, float]
            RGBA values in the range [0, 1].
        """
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
        Frame an object or MolecularEntity in the active camera view.

        Parameters
        ----------
        obj : bpy.types.Object | MolecularEntity
            Blender object or Molecular Nodes entity to frame.
        viewpoint : Viewpoints | str, optional
            Viewing direction along a principal axis. One of
            {"default", "front", "back", "top", "bottom", "left", "right"}.

        """
        if isinstance(obj, MolecularEntity):
            obj = obj.object
        # set the camera viewpoint if specified
        if viewpoint is not None:
            self.camera.set_viewpoint(viewpoint)
        # set the camera to look at the object
        blender_utils.look_at_object(obj)

    def frame_view(
        self, view: list[tuple] | MolecularEntity, viewpoint: Viewpoints | None = None
    ) -> None:
        """
        Frame one or more views of Molecular entities.
        Multiple views can be combined using ``+`` before passing the result.

        Parameters
        ----------
        view : list[tuple] | MolecularEntity
            A bounding box represented by 8 three-dimensional vertices
            ``[(x, y, z), ...]`` or an entity from which a view is derived.
        viewpoint : Viewpoints | str, optional
            Viewing direction along a principal axis. One of
            {"default", "front", "back", "top", "bottom", "left", "right"}.

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
        Clear all Molecular Nodes entities from the scene.

        Notes
        -----
        This does not modify lighting, world, or render settings.
        """
        session = get_session()
        # Iterate over a list copy to avoid modifying the dict during iteration
        for entity in list(session.entities.values()):
            session.remove(entity.uuid)

    def scene_reset(
        self,
        template: Path | str | None = "Molecular Nodes",
        engine: Cycles | EEVEE | str = "EEVEE",
    ) -> None:
        """
        Reset the scene from a template or startup file.

        Parameters
        ----------
        template : pathlib.Path | str | None, default "Molecular Nodes"
            Name of an installed Blender app template, a path to a ``.blend``
            file, or ``None`` to use Blender's default startup file.
        engine : EEVEE | Cycles | str, default "EEVEE"
            Render engine to configure after loading the template.

        Raises
        ------
        ValueError
            If ``template`` is not ``None``, not a valid ``.blend`` file path,
            and not a known app template name.
        """
        if template is None:
            bpy.ops.wm.read_homefile(app_template="")
        else:
            file = Path(template) if isinstance(template, str) else template
            if file.is_file() and file.suffix == ".blend":
                self.load(file)
            elif isinstance(template, str):
                if template in list_templates():
                    bpy.ops.wm.read_homefile(app_template=template)
                else:
                    raise ValueError(
                        f"Template '{template}' is not a valid .blend file or app template name."
                    )

            else:
                raise ValueError(
                    f"Template '{template}' is not a valid .blend file or app template name."
                )

        if engine:
            self.engine = engine

    def load(self, path: str | Path) -> None:
        """
        Load a .blend file replacing the current scene.

        Parameters
        ----------
        path : str | Path
            The file path to the .blend file to load.
        """
        file_path = Path(path) if isinstance(path, str) else path
        if not file_path.is_file() or file_path.suffix != ".blend":
            raise ValueError(f"File '{path}' is not a valid .blend file.")
        bpy.ops.wm.open_mainfile(filepath=str(file_path.resolve()))

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
        # temporary properties to override
        override_props = [
            (render_settings, "use_file_extension", True),
            (scene, "frame_current", render_frame),
        ]
        if bpy.app.version >= (5, 0, 0):
            override_props.append((image_settings, "media_type", "IMAGE"))
        override_props.append((image_settings, "file_format", file_format))
        with ExitStack() as stack:
            # set the use_file_extension to auto generate file extension
            # set the file_format to the specified one
            temp_override_properties(stack, override_props)
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
        # temporary properties to override
        override_props = [
            (render_settings, "use_lock_interface", True),
            (render_settings, "resolution_percentage", render_scale),
            (render_settings, "filepath", ""),
            (scene, "frame_current", start),
        ]
        if bpy.app.version >= (5, 0, 0):
            override_props.append((image_settings, "media_type", "IMAGE"))
        override_props.append((image_settings, "file_format", "PNG"))
        # create a temporary directory
        with tempfile.TemporaryDirectory() as tmp_dir:
            # render individual frames
            with ExitStack() as stack:
                temp_override_properties(stack, override_props)
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
            # temporary properties to override
            override_props = [
                (render_settings, "use_lock_interface", True),
                (render_settings, "resolution_percentage", render_scale),
                (render_settings, "filepath", ""),
                (scene, "frame_current", start),
                (scene, "frame_start", start),
                (scene, "frame_end", end),
            ]
            if bpy.app.version >= (5, 0, 0):
                override_props.append((image_settings, "media_type", "VIDEO"))
            override_props.append((image_settings, "file_format", "FFMPEG"))
            override_props.append((render_settings.ffmpeg, "format", "MPEG4"))
            with ExitStack() as stack:
                temp_override_properties(stack, override_props)
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
