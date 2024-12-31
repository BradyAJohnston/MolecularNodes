import bpy


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
    def resoution(self, value: tuple[int, int]) -> None:
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
            value = "BLENDER_EEVEE"
        AVAILABLE_ENGINES = ("BLENDER_EEVEE", "CYCLES", "BLENDER_WORKBENCH")
        if value in AVAILABLE_ENGINES:
            self.scene.render.engine = value  # type: ignore
        else:
            raise ValueError(
                f"Invalid render engine. Choose from: {AVAILABLE_ENGINES=}."
            )

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
