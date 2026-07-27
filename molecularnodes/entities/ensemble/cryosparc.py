from pathlib import Path
import databpy as db
import numpy as np
from .base import Ensemble, EntityType


class Dataset:
    def __init__(self, dset: np.typing.NDArray):
        self.dset = dset

    def __getitem__(self, key: str) -> np.typing.NDArray:
        return self.dset[key]

    def __contains__(self, key: str) -> bool:
        if self.dset.dtype.names is None:
            return False
        return key in self.dset.dtype.names

    def __len__(self) -> int:
        return len(self.dset)

    def get(
        self, key: str, fallback: np.typing.NDArray | None = None
    ) -> np.typing.NDArray | None:
        if key in self:
            return self[key]
        else:
            return fallback

    @property
    def alignment_psize(self) -> np.ndarray[tuple[int], np.dtype[np.float32]] | None:
        """
        Pixel size used during alignment.

        The pixel size used during alignment in CryoSPARC may differ from extracted pixel
        size. Shifts are always reported in aligned pixel size.
        """
        return self.get("alignments3D/psize_A", self.get("alignments2D/psize_A"))

    @property
    def shifts(self) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
        """
        Aligned particle shifts in A.
        """
        shifts = self.get("alignments3D/shift", self.get("alignments2D/shift"))
        if shifts is None or self.alignment_psize is None:
            return np.zeros((len(self), 2), np.float32)
        else:
            return shifts * self.alignment_psize

    @property
    def positions(self) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
        """
        Particle positions in micrograph fractions (X, Y).
        """
        xpos = self.get("location/center_x_frac")
        ypos = self.get("location/center_y_frac")
        if xpos is None or ypos is None:
            return np.zeros((len(self), 2), np.float32)
        return np.column_stack([xpos, ypos])

    @property
    def defocus(self) -> np.ndarray[tuple[int], np.dtype[np.float32]]:
        """
        Defocus in microns.
        """
        df1 = self.get("ctf/df1_A")
        df2 = self.get("ctf/df2_A")
        if df1 is None or df2 is None:
            return np.zeros(len(self), np.float32)
        return np.mean(np.column_stack([df1, df2]) / 10_000, axis=1)


class CryoSPARC(Ensemble):
    def __init__(self, file_path: Path):
        if file_path.suffix != ".cs":
            raise ValueError(f"{file_path.name} is not a CryoSPARC .cs file")
        super().__init__(file_path=file_path)
        self._entity_type = EntityType.ENSEMBLE_CRYOSPARC
        self.dset = Dataset(np.load(self.file_path))
        if len(self.dset) == 0:
            raise ValueError("Dataset has no rows")

    def create_object(
        self,
        name: str = "CryoSPARC Ensemble",
        node_setup: bool = True,
        world_scale: float = 0.01,
        fraction: float = 1.0,
        simplify=False,
    ) -> db.BlenderObject:
        defocus = self.dset.defocus.reshape((len(self.dset), -1))
        positions = np.append(self.dset.positions, defocus, axis=1)
        bob = db.create_bob(positions, name=name)
        return bob
