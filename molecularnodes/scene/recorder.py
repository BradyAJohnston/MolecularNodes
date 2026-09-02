import os
import shutil
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Literal
import bpy
from ..utils import suppress_stdout

if TYPE_CHECKING:
    from .base import Canvas

# formats an animation can be assembled into; matching is case-insensitive
_ANIMATION_FORMATS = Literal["MP4", "GIF", "mp4", "gif"]

try:
    from IPython.display import Image, Video
except ImportError:
    Image = None
    Video = None


def _png_size(path: Path) -> tuple[int, int]:
    """Read a PNG's pixel size from its IHDR chunk, which directly follows
    the 8-byte signature: 4 bytes length, 4 bytes "IHDR", then width and
    height as big-endian 32-bit integers."""
    with open(path, "rb") as f:
        header = f.read(24)
    if len(header) < 24 or header[12:16] != b"IHDR":
        raise ValueError(f"'{path}' is not a valid PNG file")
    return (
        int.from_bytes(header[16:20], "big"),
        int.from_bytes(header[20:24], "big"),
    )


def _resolve_format(
    path: str | Path | None, format: _ANIMATION_FORMATS | None
) -> Literal["MP4", "GIF"]:
    """Resolve an animation output format, inferring it from the path suffix
    when not given explicitly and validating it either way."""
    if format is None:
        is_gif = path is not None and Path(path).suffix.lower() == ".gif"
        return "GIF" if is_gif else "MP4"
    resolved = format.upper()
    if resolved not in ("MP4", "GIF"):
        raise ValueError(f"Format '{format}' must be either 'MP4' or 'GIF'")
    return resolved


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


def _encode_mp4(
    frames: list[Path], out_dir: Path, fps: float, resolution: tuple[int, int]
) -> Path:
    """Encode already-rendered frames into an MP4 without re-rendering them.

    The frames are loaded as an image strip in the sequencer of a throwaway
    scene, which Blender's bundled FFMPEG then encodes with the same container
    and codec settings that ``Canvas.animation`` uses.
    """
    scene = bpy.data.scenes.new("MN_encode_frames")
    assert scene
    try:
        render = scene.render
        view_settings = scene.view_settings
        assert render.ffmpeg and render.image_settings and view_settings
        render.resolution_x, render.resolution_y = resolution
        render.resolution_percentage = 100
        render.use_compositing = False
        render.use_sequencer = True
        # fps is an int in Blender; fps_base carries the fractional remainder
        render.fps = max(1, round(fps))
        render.fps_base = render.fps / fps
        scene.frame_start = 1
        scene.frame_end = len(frames)
        # the PNGs already carry the canvas's view transform baked in, so the
        # encode pass must not tone-map them a second time
        view_settings.view_transform = "Standard"  # ty: ignore[invalid-assignment]
        render.image_settings.media_type = "VIDEO"
        render.image_settings.file_format = "FFMPEG"
        render.ffmpeg.format = "MPEG4"
        render.ffmpeg.codec = "H264"
        # Blender appends the frame range and extension to the directory path
        render.filepath = os.path.join(str(out_dir), "")

        editor = scene.sequence_editor_create()
        assert editor
        strip = editor.strips.new_image(
            name="frames",
            filepath=str(frames[0]),
            channel=1,
            frame_start=1,
            fit_method="ORIGINAL",
        )
        assert isinstance(strip, bpy.types.ImageStrip)
        # each appended element extends the strip by one frame
        for frame in frames[1:]:
            strip.elements.append(frame.name)

        with bpy.context.temp_override(scene=scene):
            with suppress_stdout():
                bpy.ops.render.render(animation=True)
    finally:
        bpy.data.scenes.remove(scene)
    return next(out_dir.glob("*.mp4"))


class FrameRecorder:
    """
    Collect manually rendered frames into an animation.

    Where [](`~mn.Canvas.animation`) renders an animation by playing back the
    Blender timeline, a recorder gives the loop to you: change anything about
    the scene between frames - the timeline, the camera, styles, colors - and
    call :meth:`render` to capture the scene as it stands as the next frame.
    :meth:`finalize` assembles the captured frames into an MP4 or GIF.

    Created via [](`~mn.Canvas.record`) rather than directly. Frames are
    rendered as PNGs into a temporary directory that lives (and is cleaned up)
    with the recorder, or into ``frames_dir`` when one is given, where they
    stay after the recorder is gone. With ``overwrite=False``, frames already
    in ``frames_dir`` are reused rather than re-rendered, resuming an
    interrupted recording.

    A recorder is also a context manager: when a ``path`` was given to
    [](`~mn.Canvas.record`), leaving the ``with`` block without an exception
    finalizes to that path automatically.

    Examples
    --------
    Explicit finalize, which returns the animation for notebook display::

        movie = canvas.record(fps=24)
        for i in range(120):
            traj.frame = i
            canvas.look_at(traj)
            movie.render()
        movie.finalize("wobble.mp4")

    As a context manager, finalizing on exit::

        with canvas.record("wobble.mp4", fps=24) as movie:
            for i in range(120):
                traj.frame = i
                movie.render()

    See Also
    --------
    molecularnodes.Canvas.record : Create a recorder bound to a canvas.
    molecularnodes.Canvas.animation : Render an animation from the timeline.
    """

    def __init__(
        self,
        canvas: "Canvas",
        path: str | Path | None = None,
        fps: float | None = None,
        render_scale: int = 100,
        frames_dir: str | Path | None = None,
        overwrite: bool = True,
    ) -> None:
        self._canvas = canvas
        self._path = Path(path) if path is not None else None
        self._fps = fps
        self._render_scale = render_scale
        self._overwrite = overwrite
        if frames_dir is None:
            if not overwrite:
                raise ValueError(
                    "overwrite=False requires frames_dir - without a "
                    "persistent frames directory there are no earlier frames "
                    "to resume from."
                )
            self._tmp_dir = tempfile.TemporaryDirectory(prefix="mn_frames_")
            self._frames_dir = Path(self._tmp_dir.name)
        else:
            # frames written to a user-given directory outlive the recorder
            self._tmp_dir = None
            self._frames_dir = Path(frames_dir)
            self._frames_dir.mkdir(parents=True, exist_ok=True)
        self._frames: list[Path] = []
        self._pixel_size: tuple[int, int] | None = None
        self._finalized = False

    @property
    def frames(self) -> list[Path]:
        """
        Get the frames captured so far.

        Returns
        -------
        list[pathlib.Path]
            Paths of the rendered PNG frames, in capture order. Without a
            ``frames_dir``, the files live in a temporary directory that is
            removed with the recorder, so copy them elsewhere to keep them
            beyond it.
        """
        return list(self._frames)

    def __len__(self) -> int:
        return len(self._frames)

    def __repr__(self) -> str:
        size = "" if self._pixel_size is None else f" at {self._pixel_size}"
        return f"<FrameRecorder: {len(self._frames)} frames{size}>"

    def _frame_size(self, render_scale: int) -> tuple[int, int]:
        # Blender truncates when applying the resolution percentage
        x, y = self._canvas.resolution
        return (x * render_scale // 100, y * render_scale // 100)

    def render(self, render_scale: int | None = None) -> Path:
        """
        Render the scene as it stands and store it as the next frame.

        Parameters
        ----------
        render_scale : int, optional
            Scale of the rendered frame with respect to the resolution,
            overriding the recorder's default for this frame.

        Returns
        -------
        pathlib.Path
            Path of the rendered PNG frame.

        Raises
        ------
        ValueError
            If the frame's pixel size differs from the frames already
            captured - all frames of an animation must share one resolution.

        Notes
        -----
        With ``overwrite=False`` (and a ``frames_dir``), a frame whose
        numbered file already exists is reused instead of re-rendered, so
        re-running the same loop after a crash only renders the frames that
        are missing.
        """
        scale = self._render_scale if render_scale is None else render_scale
        stem = self._frames_dir / f"{len(self._frames):05d}"
        existing = stem.with_name(stem.name + ".png")
        if not self._overwrite and existing.exists():
            # resume: take the frame's size from the file itself, so frames
            # rendered by an earlier run are checked as they actually are
            frame, size = existing, _png_size(existing)
        else:
            frame, size = None, self._frame_size(scale)
        if self._pixel_size is None:
            self._pixel_size = size
        elif size != self._pixel_size:
            raise ValueError(
                f"Frame size {size} does not match the recording's "
                f"{self._pixel_size}; all frames of an animation must share "
                "one resolution."
            )
        if frame is None:
            frame = self._canvas._render_still(
                stem, file_format="PNG", render_scale=scale
            )
        self._frames.append(frame)
        return frame

    def finalize(
        self,
        path: str | Path | None = None,
        fps: float | None = None,
        format: _ANIMATION_FORMATS | None = None,
    ) -> "Video | Image | None":
        """
        Assemble the captured frames into an animation.

        The recorder is left intact: more frames can be rendered and
        ``finalize`` called again, e.g. to write both an MP4 and a GIF of the
        same recording.

        Parameters
        ----------
        path : str | Path | None, optional
            File path to write the animation to, defaulting to the one given
            to [](`~mn.Canvas.record`). The animation is returned for display
            regardless of whether a path is given.
        fps : float, optional
            Frame rate of the animation. Defaults to the rate given to
            [](`~mn.Canvas.record`), or failing that the scene's fps.
        format : str, optional
            Output format, either ``"MP4"`` or ``"GIF"`` (case-insensitive).
            When not specified, inferred from the suffix of ``path``,
            defaulting to MP4. GIF output requires the ``pillow`` package.

        Returns
        -------
        IPython.display.Video | IPython.display.Image | None
            The assembled animation (a ``Video`` for MP4, an ``Image`` for
            GIF), which displays automatically as the result of a notebook
            cell. ``None`` if IPython is not installed.

        Raises
        ------
        RuntimeError
            If no frames have been rendered yet.
        """
        if not self._frames:
            raise RuntimeError(
                "No frames have been rendered - call render() before finalize()."
            )
        if path is None:
            path = self._path
        if fps is None:
            fps = self._fps if self._fps is not None else self._canvas.fps
        format = _resolve_format(path, format)
        assert self._pixel_size

        with tempfile.TemporaryDirectory() as tmp_dir:
            if format == "MP4":
                output_file = _encode_mp4(
                    self._frames, Path(tmp_dir), fps, self._pixel_size
                )
            else:
                output_file = Path(tmp_dir) / "animation.gif"
                _write_gif(self._frames, output_file, fps)
            data = output_file.read_bytes()
            if path:
                shutil.copy(output_file, path)
        self._finalized = True

        if format == "MP4":
            if Video is None:
                return None
            return Video(data=data, embed=True, mimetype="video/mp4")
        if Image is None:
            return None
        return Image(data=data, format="gif")

    def __enter__(self) -> "FrameRecorder":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        # only a clean exit finalizes, and only when a destination was given
        # up front and finalize() hasn't already been called in the block
        if exc_type is None and self._path is not None and not self._finalized:
            self.finalize()
