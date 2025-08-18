import os
import shutil
import tempfile
from contextlib import ExitStack
from pathlib import Path
import bpy
from tqdm.auto import tqdm
from ..entities.base import MolecularEntity
from ..utils import suppress_stdout, temp_override_properties
from .engines import EEVEE, Cycles

try:
    from IPython.display import Image, Video, display
except ImportError:
    Image = None
    Video = None
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
            "BLENDER_EEVEE",
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

    def frame_object(self, obj: bpy.types.Object | MolecularEntity) -> None:
        if isinstance(obj, MolecularEntity):
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
        scene = bpy.context.scene
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

        scene = bpy.context.scene
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
