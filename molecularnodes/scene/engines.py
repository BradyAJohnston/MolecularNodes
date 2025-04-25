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
    def __init__(self, samples: int = 64):
        self._name = "BLENDER_EEVEE_NEXT"
        self.samples = samples
        self._enable_engine()

    @property
    def engine(self):
        return bpy.context.scene.eevee

    @property
    def samples(self) -> int:
        return self.engine.taa_render_samples

    @samples.setter
    def samples(self, value: int) -> None:
        self.engine.taa_render_samples = value


class Cycles(RenderEngine):
    def __init__(self, samples: int = 256, device: str = "GPU"):
        self._name = "CYCLES"
        self.samples = samples
        self.device = device
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
            enable_optimal_gpu()
