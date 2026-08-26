import warnings
from abc import ABC
import bpy
from .render import enable_optimal_gpu


class RenderEngine(ABC):
    """
    Configuration for a render engine.

    Constructing an engine stores the requested settings without touching the
    scene; they are applied when the engine is activated by assigning it to
    ``Canvas.engine``. The properties always read from and write to the
    active scene directly.
    """

    _name = "RenderEngineName"

    @property
    def name(self) -> str:
        return self._name

    def _enable_engine(self):
        bpy.context.scene.render.engine = self.name  # type: ignore

    def _apply_settings(self) -> None:
        """Write the settings given at construction to the scene."""


class EEVEE(RenderEngine):
    _name = "BLENDER_EEVEE"

    def __init__(self, samples: int = 64, raytracing: bool = True):
        self._init_samples = samples
        self._init_raytracing = raytracing

    def _apply_settings(self) -> None:
        self.samples = self._init_samples
        self.raytracing = self._init_raytracing

    @property
    def engine(self):
        return bpy.context.scene.eevee

    @property
    def samples(self) -> int:
        return self.engine.taa_render_samples

    @samples.setter
    def samples(self, value: int) -> None:
        self.engine.taa_render_samples = value

    @property
    def raytracing(self) -> bool:
        return self.engine.use_raytracing

    @raytracing.setter
    def raytracing(self, value: bool) -> None:
        self.engine.use_raytracing = value


class Cycles(RenderEngine):
    _name = "CYCLES"

    def __init__(
        self,
        samples: int = 256,
        device: str = "GPU",
        denoise: bool = True,
        denoise_gpu: bool = True,
    ):
        self._init_samples = samples
        self._init_device = device
        self._init_denoise = denoise
        self._init_denoise_gpu = denoise_gpu

    def _apply_settings(self) -> None:
        self.samples = self._init_samples
        self.device = self._init_device
        self.denoise = self._init_denoise
        self.denoise_gpu = self._init_denoise_gpu

    @property
    def engine(self):
        return bpy.context.scene.cycles

    @property
    def samples(self) -> int:
        return self.engine.samples

    @samples.setter
    def samples(self, value: int) -> None:
        self.engine.samples = value

    @property
    def device(self) -> str:
        return self.engine.device

    @device.setter
    def device(self, value: str):
        value = value.upper()
        self.engine.device = value
        if value == "GPU":
            try:
                enable_optimal_gpu()
            except RuntimeError:
                warnings.warn(
                    "Failed to enable GPU, defaulting back to CPU render device"
                )
                self.engine.device = "CPU"

    @property
    def denoise(self) -> bool:
        return self.engine.use_denoising

    @denoise.setter
    def denoise(self, value: bool) -> None:
        self.engine.use_denoising = value

    @property
    def denoise_gpu(self) -> bool:
        return self.engine.denoising_use_gpu

    @denoise_gpu.setter
    def denoise_gpu(self, value: bool) -> None:
        self.engine.denoising_use_gpu = value
