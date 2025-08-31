import warnings
from abc import ABC
import bpy
from .render import enable_optimal_gpu


class RenderEngine(ABC):
    _name = "RenderEngineName"

    @property
    def name(self) -> str:
        return self._name

    def _enable_engine(self):
        bpy.context.scene.render.engine = self.name  # type: ignore


class EEVEE(RenderEngine):
    def __init__(self, samples: int = 64, raytracing: bool = True):
        # TODO: remove this check when Blender 5.0 is the minimum version
        if bpy.app.version >= (5, 0, 0):
            self._name = "BLENDER_EEVEE"
        else:
            self._name = "BLENDER_EEVEE_NEXT"
        self.samples = samples
        self.raytracing = raytracing

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
    def __init__(
        self,
        samples: int = 256,
        device: str = "GPU",
        denoise: bool = True,
        denoise_gpu: bool = True,
    ):
        self._name = "CYCLES"
        self.samples = samples
        self.device = device
        self.denoise = denoise
        self.denoise_gpu = denoise_gpu
        self._enable_engine()

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
            except TypeError:
                warnings.warn(
                    "Failed to enable GPU, defaulting back to CPU render device"
                )
                self.engine.device = "CPU"

    @property
    def denoise(self) -> bool:
        return self.engine.use_denoising

    @denoise.setter
    def denoise(self, value: bool) -> None:
        self.engine.use_denoiseing = value

    @property
    def denoise_gpu(self) -> bool:
        return self.engine.denoising_use_gpu

    @denoise_gpu.setter
    def denoise_gpu(self, value: bool) -> None:
        self.engine.denoising_use_gpu = value
