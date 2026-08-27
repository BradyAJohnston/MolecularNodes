import os
import shutil
import tempfile
from contextlib import ExitStack, contextmanager
from enum import StrEnum
from pathlib import Path
from typing import Literal, Sequence, Tuple
import bpy
from tqdm.auto import tqdm
from .. import assets
from ..assets.template import list_templates
from ..blender import utils as blender_utils
from ..entities.base import MolecularEntity
from ..session import get_session
from ..ui import addon
from ..utils import suppress_stdout, temp_override_properties
from .camera import Camera, Viewpoint
from .compositor import CompositorTree, setup_compositor
from .engines import EEVEE, Cycles
from .world import WorldTree

try:
    from IPython.display import Image, Video
except ImportError:
    Image = None
    Video = None


class ViewTransform(StrEnum):
    STANDARD = "Standard"
    KHRONOS = "Khronos PBR Neutral"
    AGX = "AgX"
    FILMIC = "Filmic"
    FILMIC_LOG = "Filmic Log"
    FALSE_COLOR = "False Color"
    RAW = "Raw"

    @classmethod
    def _missing_(cls, value: object) -> "ViewTransform | None":
        # called only when the exact-value lookup fails; allow the full name or the
        # short enum name, case-insensitively (e.g. "agx", "filmic_log")
        if isinstance(value, str):
            key = value.strip().lower()
            for member in cls:
                if key in (member.value.lower(), member.name.lower()):
                    return member
        return None


_RENDER_ENGINES = Literal["EEVEE", "CYCLES"]

# render passes commonly needed for compositor effects, mapped to the
# `bpy.types.ViewLayer` toggle that enables each one
RENDER_PASSES = {
    "combined": "use_pass_combined",
    "z": "use_pass_z",
    "mist": "use_pass_mist",
    "normal": "use_pass_normal",
    "position": "use_pass_position",
    "vector": "use_pass_vector",
    "diffuse_color": "use_pass_diffuse_color",
    "emit": "use_pass_emit",
    "environment": "use_pass_environment",
    "ambient_occlusion": "use_pass_ambient_occlusion",
}

# still-image formats that IPython can embed inline, mapped to the format
# name expected by IPython.display.Image
_IPYTHON_IMAGE_FORMATS = {"PNG": "png", "JPEG": "jpeg"}


@contextmanager
def _restore_frame(scene: bpy.types.Scene):
    """Restore the scene's current frame (via ``frame_set``, so animation data
    is re-evaluated) on exit."""
    original = scene.frame_current
    try:
        yield
    finally:
        scene.frame_set(original)


@contextmanager
def _render_progress(total: int, desc: str):
    """Show a tqdm progress bar driven by the ``render_write`` handler, which
    fires once per frame written during an animation render."""
    bar = tqdm(total=total, desc=desc)

    def _on_frame_write(*args) -> None:
        bar.update(1)

    bpy.app.handlers.render_write.append(_on_frame_write)
    try:
        yield
    finally:
        bpy.app.handlers.render_write.remove(_on_frame_write)
        bar.close()


def _write_gif(frames: list[Path], path: Path, fps: float) -> None:
    """Assemble rendered frames into an animated GIF using Pillow."""
    try:
        from PIL import Image as PILImage
    except ImportError as e:
        raise ImportError(
            "GIF output requires the optional dependency 'pillow'. "
            "Install it with `pip install pillow`."
        ) from e
    images = [PILImage.open(frame) for frame in frames]
    images[0].save(
        path,
        save_all=True,
        append_images=images[1:],
        duration=round(1000 / fps),
        loop=0,
        disposal=2,
    )


class Canvas:
    """
    High-level render controller for Blender scenes.

    Canvas configures the active Blender scene for Molecular Nodes renders
    (engine, resolution, transparency, color management), exposes convenient
    properties for common render settings, and provides helpers to point the
    camera at objects/views and render stills or animations.

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
    compositor : molecularnodes.scene.compositor.CompositorTree
        Builder for the scene compositor node tree (post-processing effects).
    world : molecularnodes.scene.world.WorldTree
        Builder for the world shader node tree (lighting and background).
    samples : int
        Render sample count on the active engine.
    frame : int
        Current scene frame.
    frame_range : tuple[int, int]
        Scene ``(frame_start, frame_end)``.
    render_scale : int
        Render resolution percentage (100 = full).
    exposure : float
        Color-management exposure.
    gamma : float
        Color-management gamma.
    look : str
        Color-management look (contrast preset).
    passes : list[str]
        Enabled render passes (see :data:`RENDER_PASSES`).

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
        engine: EEVEE | Cycles | _RENDER_ENGINES = "EEVEE",
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
        self._compositor: CompositorTree | None = None
        self._world: WorldTree | None = None
        if template:
            self.scene_reset(template=template)
        else:
            self._scene_changed()
        self.engine = engine
        self.resolution = resolution
        self.camera = Camera()
        self.transparent = transparent

    def _scene_changed(self) -> None:
        """
        Rebind scene-backed state after the scene has been replaced.

        The world and compositor builders hold references to the old scene's
        node trees, so they are discarded and the compositor (with the
        annotation overlay) is prepared on the new scene.
        """
        self._compositor = None
        self._world = None
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
        scene = bpy.context.scene
        assert scene
        return scene

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
        Get the active render engine.

        Returns
        -------
        EEVEE | Cycles
            The active engine configuration object. Kept in sync with the
            scene, so an engine switch made outside of this class (e.g. in the
            Blender UI or by loading a .blend file) is reflected here.
        """
        engine_id = self.scene.render.engine
        if self._engine.name != engine_id:
            # the engine was changed outside of this class - rebind without
            # applying the constructor defaults over the scene's settings
            self._engine = Cycles() if engine_id == "CYCLES" else EEVEE()
        return self._engine

    @engine.setter
    def engine(self, value: Cycles | EEVEE | _RENDER_ENGINES) -> None:
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
        # settings given to the engine's constructor are applied on activation
        self._engine._apply_settings()

    @property
    def fps(self) -> int:
        """
        Get and set the frame rate expressed in frames per second.

        Returns
        -------
        int
            The current FPS value.
        """
        return self.scene.render.fps

    @fps.setter
    def fps(self, value: int) -> None:
        """
        Get and set the frame rate expressed in frames per second.

        Parameters
        ----------
        value : int
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
    def compositor(self) -> CompositorTree:
        """
        The scene compositor node tree.

        Build post-processing effects with ``nodebpy.compositor`` nodes, either
        by appending (``with canvas.compositor as tree``) or from a clean graph
        (``canvas.compositor.reset()``).

        Returns
        -------
        molecularnodes.scene.compositor.CompositorTree
            Builder bound to the active scene's compositor tree.
        """
        if self._compositor is None:
            self._compositor = CompositorTree(self.scene)
        return self._compositor

    @property
    def world(self) -> WorldTree:
        """
        The scene world shader node tree (lighting & background).

        Returns
        -------
        molecularnodes.scene.world.WorldTree
            Builder bound to the active scene's world shader tree.
        """
        if self._world is None:
            self._world = WorldTree(self.scene)
        return self._world

    @property
    def background(self) -> Tuple[float, float, float, float]:
        """
        Get the world background color.

        Returns
        -------
        tuple[float, float, float, float]
            RGBA values in the range [0, 1].
        """
        return self.world.background

    @background.setter
    def background(self, value: Tuple[float, float, float, float]) -> None:
        """
        Set the world background color.

        Parameters
        ----------
        value : tuple[float, float, float, float]
            RGBA values in the range [0, 1].
        """
        self.world.background = value

    @property
    def samples(self) -> int:
        """
        Get the render sample count from the active engine.

        Returns
        -------
        int
            Number of samples the active render engine is configured for.
        """
        return self.engine.samples

    @samples.setter
    def samples(self, value: int) -> None:
        """
        Set the render sample count on the active engine.

        Parameters
        ----------
        value : int
            Number of samples to render with.
        """
        self.engine.samples = value

    @property
    def frame(self) -> int:
        """
        Get the current scene frame.

        Returns
        -------
        int
            The current frame number.
        """
        return self.scene.frame_current

    @frame.setter
    def frame(self, value: int) -> None:
        """
        Set the current scene frame.

        Parameters
        ----------
        value : int
            The frame number to set. Uses ``frame_set`` so animation data updates.
        """
        self.scene.frame_set(value)

    @property
    def frame_range(self) -> tuple[int, int]:
        """
        Get the scene frame range.

        Returns
        -------
        tuple[int, int]
            The ``(frame_start, frame_end)`` of the scene.
        """
        return (self.scene.frame_start, self.scene.frame_end)

    @frame_range.setter
    def frame_range(self, value: tuple[int, int]) -> None:
        """
        Set the scene frame range.

        Parameters
        ----------
        value : tuple[int, int]
            The ``(frame_start, frame_end)`` to set.
        """
        self.scene.frame_start, self.scene.frame_end = value

    @property
    def render_scale(self) -> int:
        """
        Get the render resolution percentage.

        Returns
        -------
        int
            The resolution percentage applied to the render (100 = full).
        """
        return self.scene.render.resolution_percentage

    @render_scale.setter
    def render_scale(self, value: int) -> None:
        """
        Set the render resolution percentage.

        Parameters
        ----------
        value : int
            Percentage of the resolution to render at (100 = full).
        """
        self.scene.render.resolution_percentage = value

    @property
    def _view_settings(self) -> bpy.types.ColorManagedViewSettings:

        view_settings = self.scene.view_settings
        assert view_settings
        return view_settings

    @property
    def exposure(self) -> float:
        """
        Get the color-management exposure.

        Returns
        -------
        float
            Exposure applied in the view transform.
        """
        return self._view_settings.exposure

    @exposure.setter
    def exposure(self, value: float) -> None:
        """
        Set the color-management exposure.

        Parameters
        ----------
        value : float
            Exposure to apply in the view transform.
        """
        self._view_settings.exposure = value

    @property
    def gamma(self) -> float:
        """
        Get the color-management gamma.

        Returns
        -------
        float
            Gamma applied in the view transform.
        """
        return self._view_settings.gamma

    @gamma.setter
    def gamma(self, value: float) -> None:
        """
        Set the color-management gamma.

        Parameters
        ----------
        value : float
            Gamma to apply in the view transform.
        """
        self._view_settings.gamma = value

    @property
    def look(self) -> str:
        """
        Get the color-management look (contrast preset).

        Returns
        -------
        str
            The active look, e.g. ``"None"`` or ``"AgX - Medium High Contrast"``.
        """
        return self._view_settings.look

    @look.setter
    def look(self, value: str) -> None:
        """
        Set the color-management look (contrast preset).

        Parameters
        ----------
        value : str
            The look to apply. Valid values depend on the active view transform.
        """
        self._view_settings.look = value  # ty: ignore[invalid-assignment]

    @property
    def passes(self) -> list[str]:
        """
        Get the enabled render passes.

        Returns
        -------
        list[str]
            Names of the enabled passes from :data:`RENDER_PASSES`.
        """
        view_layer = self.scene.view_layers[0]
        return [
            name for name, attr in RENDER_PASSES.items() if getattr(view_layer, attr)
        ]

    @passes.setter
    def passes(self, value: list[str]) -> None:
        """
        Enable a set of render passes (for use in the compositor).

        Parameters
        ----------
        value : list[str]
            Names of passes to enable from :data:`RENDER_PASSES`; all others are
            disabled.

        Raises
        ------
        ValueError
            If a name is not a recognised pass.
        """
        view_layer = self.scene.view_layers[0]
        requested = set(value)
        unknown = requested - set(RENDER_PASSES)
        if unknown:
            raise ValueError(
                f"Unknown render pass(es): {sorted(unknown)}. "
                f"Valid passes are {sorted(RENDER_PASSES)}."
            )
        for name, attr in RENDER_PASSES.items():
            setattr(view_layer, attr, name in requested)

    @property
    def view_transform(self) -> ViewTransform:
        """
        Get the current view transform setting.

        Returns
        -------
        ViewTransform
            The current view transform. As a ``StrEnum`` this also compares equal
            to its Blender string value.
        """
        return ViewTransform(self._view_settings.view_transform)

    @view_transform.setter
    def view_transform(self, value: ViewTransform) -> None:
        """
        Set the view transform setting.

        Parameters
        ----------
        value : ViewTransform | str
            The view transform to set. Accepts a ``ViewTransform``, its full name
            (case-insensitive), or its short enum name (e.g. ``"agx"``,
            ``"filmic_log"``).
        """
        self._view_settings.view_transform = ViewTransform(value)  # ty: ignore[invalid-assignment]

    def look_at(
        self,
        target: MolecularEntity | bpy.types.Object | list[tuple],
        viewpoint: Viewpoint | str | Sequence[float] | None = None,
    ) -> None:
        """
        Position the camera to look at and contain a target.

        Parameters
        ----------
        target : MolecularEntity | bpy.types.Object | list[tuple]
            What to look at: a Molecular Nodes entity (via its current view),
            a Blender object, or a bounding box of 8 three-dimensional
            vertices ``[(x, y, z), ...]`` as returned by ``get_view()``.
            Multiple views can be combined with ``+`` before passing.
        viewpoint : Viewpoint | str | Sequence[float], optional
            Viewing direction along a principal axis — one of
            {"default", "front", "back", "top", "bottom", "left", "right"} —
            or a custom XYZ Euler rotation as three angles in degrees.

        """
        # set the camera viewpoint if specified
        if viewpoint is not None:
            self.camera.set_viewpoint(viewpoint)
        if isinstance(target, MolecularEntity):
            target = target.get_view()
        if isinstance(target, bpy.types.Object):
            blender_utils.look_at_object(target)
        else:
            blender_utils.look_at_bbox(target)

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
        engine: Cycles | EEVEE | _RENDER_ENGINES = "EEVEE",
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
            self._scene_changed()
        else:
            file = Path(template) if isinstance(template, str) else template
            if file.is_file() and file.suffix == ".blend":
                self.load(file)
            elif isinstance(template, str):
                if template in list_templates():
                    bpy.ops.wm.read_homefile(app_template=template)
                    self._scene_changed()
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
        self._scene_changed()

    def snapshot(
        self,
        path: str | Path | None = None,
        frame: int | None = None,
        file_format: str = "PNG",
        render_scale: int = 100,
    ) -> "Image | None":
        """
        Render an image of the current scene.

        Parameters
        ----------
        path : str | Path | None, optional
            File path to write the rendered image to. The image is returned
            for display regardless of whether a path is given.

        frame : int, optional
            Frame number of scene to render. When not specified,
            current scene's current_frame is used

        file_format : str, optional
            File format of the rendered image.

        render_scale : int, optional
            Scale of the rendered image with respect to the resolution.

        Returns
        -------
        IPython.display.Image | None
            The rendered image, which displays automatically as the result of
            a notebook cell. ``None`` if IPython is not installed or the
            format cannot be displayed in a notebook (e.g. ``"OPEN_EXR"``).
        """
        scene = self.scene
        render_settings = scene.render
        image_settings = render_settings.image_settings
        # temporary properties to override
        override_props = [
            (image_settings, "media_type", "IMAGE"),
            (image_settings, "file_format", file_format),
            (render_settings, "resolution_percentage", render_scale),
        ]
        with ExitStack() as stack:
            temp_override_properties(stack, override_props)
            # restore the current frame (and animation state) afterwards
            stack.enter_context(_restore_frame(scene))
            if frame is not None:
                # only frame_set will update animation data,
                # scene.frame_current will only update the timeline
                scene.frame_set(frame)
            # render into a temporary directory, with the extension matching
            # the requested file format
            tmp_dir = stack.enter_context(tempfile.TemporaryDirectory())
            render_file = Path(tmp_dir) / f"snapshot{render_settings.file_extension}"
            temp_override_properties(
                stack, [(render_settings, "filepath", str(render_file))]
            )
            # suppress stdout output generated by render process
            with suppress_stdout():
                bpy.ops.render.render(write_still=True, animation=False)
            data = render_file.read_bytes()
            if path:
                shutil.copy(render_file, path)
        ipython_format = _IPYTHON_IMAGE_FORMATS.get(file_format.upper())
        if Image is None or ipython_format is None:
            return None
        return Image(data=data, format=ipython_format)

    def animation(
        self,
        path: str | Path | None = None,
        frame_start: int | None = None,
        frame_end: int | None = None,
        render_scale: int = 100,
        fps: float | None = None,
        format: str | None = None,
    ) -> "Video | Image | None":
        """
        Render an animation of the current scene.

        Parameters
        ----------
        path : str | Path | None, optional
            File path to write the rendered animation to. The animation is
            returned for display regardless of whether a path is given.

        frame_start : int, optional
            Start frame of the animation. When not specified, current scene's
            start frame is used

        frame_end : int, optional
            End frame of the animation. When not specified, current scene's
            end frame is used

        render_scale : int, optional
            Scale of the rendered animation frames with respect to the resolution.

        fps : float, optional
            Frame rate of the animation. When not specified, the scene's
            fps is used.

        format : str, optional
            Output format, either ``"MP4"`` or ``"GIF"`` (case-insensitive).
            When not specified, inferred from the suffix of ``path``,
            defaulting to MP4. GIF output requires the ``pillow`` package.

        Returns
        -------
        IPython.display.Video | IPython.display.Image | None
            The rendered animation (a ``Video`` for MP4, an ``Image`` for
            GIF), which displays automatically as the result of a notebook
            cell. ``None`` if IPython is not installed.
        """
        # determine frame range
        start = self.frame_start if frame_start is None else frame_start
        end = self.frame_end if frame_end is None else frame_end
        if end < start:
            raise ValueError(f"End frame {end} cannot be less than start frame {start}")

        if format is None:
            is_gif = path is not None and Path(path).suffix.lower() == ".gif"
            format = "GIF" if is_gif else "MP4"
        format = format.upper()
        if format not in ("MP4", "GIF"):
            raise ValueError(f"Format '{format}' must be either 'MP4' or 'GIF'")

        scene = self.scene
        render_settings = scene.render
        image_settings = render_settings.image_settings
        # temporary properties to override
        override_props = [
            (render_settings, "use_lock_interface", True),
            (render_settings, "use_file_extension", True),
            (render_settings, "resolution_percentage", render_scale),
            (scene, "frame_start", start),
            (scene, "frame_end", end),
        ]
        if fps is not None:
            override_props.append((render_settings, "fps", fps))
        if format == "MP4":
            override_props += [
                (image_settings, "media_type", "VIDEO"),
                (image_settings, "file_format", "FFMPEG"),
                (render_settings.ffmpeg, "format", "MPEG4"),
                (render_settings.ffmpeg, "codec", "H264"),
            ]
        else:
            # GIF frames are rendered as PNGs and assembled with pillow
            override_props += [
                (image_settings, "media_type", "IMAGE"),
                (image_settings, "file_format", "PNG"),
            ]

        with tempfile.TemporaryDirectory() as tmp_dir, ExitStack() as stack:
            temp_override_properties(stack, override_props)
            # restore the current frame (and animation state) afterwards
            stack.enter_context(_restore_frame(scene))
            # Blender appends the frame range (MP4) or frame number (PNG) and
            # the extension to the output path
            temp_override_properties(
                stack, [(render_settings, "filepath", os.path.join(tmp_dir, ""))]
            )
            with _render_progress(end - start + 1, "Rendering frames"):
                with suppress_stdout():
                    bpy.ops.render.render(animation=True)

            if format == "MP4":
                output_file = next(Path(tmp_dir).glob("*.mp4"))
            else:
                # frame numbers make up the filenames, so sort numerically
                frames = sorted(Path(tmp_dir).glob("*.png"), key=lambda p: int(p.stem))
                output_file = Path(tmp_dir) / "animation.gif"
                _write_gif(frames, output_file, fps if fps is not None else self.fps)

            data = output_file.read_bytes()
            if path:
                shutil.copy(output_file, path)

        if format == "MP4":
            if Video is None:
                return None
            return Video(data=data, embed=True, mimetype="video/mp4")
        if Image is None:
            return None
        return Image(data=data, format="gif")
