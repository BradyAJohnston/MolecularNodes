from pathlib import Path
from typing import TYPE_CHECKING, overload
import databpy as db
import numpy as np
from scipy.spatial.transform import Rotation
from .base import Ensemble, EntityType

if TYPE_CHECKING:
    from typing import Any, Literal, TypeVar

    type BobField = Literal["data", "name", "atype", "domain"]
    F = TypeVar("F")


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

    @overload
    def get(self, key: str, fallback: None = None) -> np.typing.NDArray | None: ...
    @overload
    def get(self, key: str, fallback: np.typing.NDArray) -> np.typing.NDArray: ...
    def get(self, key, fallback=None):
        if key in self:
            return self[key]
        else:
            return fallback

    @property
    def psize_2d(self) -> np.ndarray[tuple[int], np.dtype[np.float32]] | None:
        """
        Pixel size used during alignment.

        The pixel size used during alignment in CryoSPARC may differ from extracted pixel
        size. Shifts are always reported in aligned pixel size.
        """
        return self.get("alignments2D/psize_A")

    @property
    def psize_3d(self) -> np.ndarray[tuple[int], np.dtype[np.float32]] | None:
        """
        Pixel size used during alignment.

        The pixel size used during alignment in CryoSPARC may differ from extracted pixel
        size. Shifts are always reported in aligned pixel size.
        """
        return self.get("alignments3D/psize_A")

    @property
    def shift2d(self) -> np.ndarray[tuple[int, int], np.dtype[np.float32]] | None:
        """
        Shifts aligned during 2D Classification
        """
        if (
            shifts := self.get("alignments2D/shift")
        ) is not None and self.psize_2d is not None:
            return shifts.astype(np.float32) * self.psize_2d.reshape(len(self), -1)

    @property
    def shift3d(self) -> np.ndarray[tuple[int, int], np.dtype[np.float32]] | None:
        """
        Aligned particle shifts in A.
        """
        if (
            shifts := self.get("alignments3D/shift")
        ) is not None and self.psize_3d is not None:
            return shifts * self.psize_3d.reshape(len(self), -1)

    @property
    def defocus(self) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
        """
        Defocus in microns.
        """
        df1 = self.get("ctf/df1_A")
        df2 = self.get("ctf/df2_A")
        if df1 is None or df2 is None:
            return np.zeros((len(self), 2), np.float32)
        return np.column_stack([df1, df2]) / 10_000

    @property
    def positions(self) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
        """
        Particle positions in micrograph fractions (X, Y).
        """
        xpos = self.get("location/center_x_frac", np.zeros(len(self), dtype=np.float32))
        ypos = self.get("location/center_y_frac", np.zeros_like(xpos))
        zpos = np.mean(self.defocus, axis=1, keepdims=True)
        return np.column_stack([xpos, ypos, zpos])

    @property
    def rotations(self) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
        """
        Pose 3D
        """
        if (poses := self.get("alignments3D/pose")) is None:
            return np.zeros((len(self), 4), np.float32)
        return Rotation.from_rotvec(poses).as_quat()


class CryoSPARC(Ensemble):
    def __init__(self, file_path: Path):
        if file_path.suffix != ".cs":
            raise ValueError(f"{file_path.name} is not a CryoSPARC .cs file")
        super().__init__(file_path=file_path)
        self._entity_type = EntityType.ENSEMBLE_CRYOSPARC
        self.dset = Dataset(np.load(self.file_path))
        if len(self.dset) == 0:
            raise ValueError("Dataset has no rows")

    def _construct_fields(self) -> "list[dict[BobField, Any]]":
        """
        Construct named attributes from CryoSPARC fields.

        Return a sorted list of funciton argument dicts.
        We sort them here so they are presented in the same order they would be in
        a CryoSPARC dataset.
        """
        fields = []
        if not np.allclose((defocus := self.dset.defocus), 0.0):
            fields.append(dict(data=defocus[:, 0], name="ctf/df1_um", atype="FLOAT"))
            fields.append(dict(data=defocus[:, 1], name="ctf/df2_um", atype="FLOAT"))
        if (shift3d := self.dset.shift3d) is not None:
            fields.append(
                dict(
                    data=shift3d,
                    name="alignments3D/shift",
                    atype="FLOAT2",
                )
            )
        if (shift2d := self.dset.shift2d) is not None:
            fields.append(
                dict(
                    data=shift2d,
                    name="alignments2D/shift",
                    atype="FLOAT2",
                )
            )

        # type checking complains about no matching overload here, but it's wrong
        return sorted(fields, key=lambda field: field.get("name"))  # type: ignore

    def create_object(
        self,
        name: str = "CryoSPARC Ensemble",
        node_setup: bool = True,
        world_scale: float = 0.01,
        fraction: float = 1.0,
        simplify=False,
    ) -> db.BlenderObject:
        # build positions and rotation outside _construct_fields because we want
        # them to appear first, like they would with a normal Blender object.
        bob = db.create_bob(self.dset.positions, name=name)
        bob.store_named_attribute(
            data=self.dset.rotations, name="rotation", atype="QUATERNION"
        )
        for field in self._construct_fields():
            bob.store_named_attribute(**field)
        if node_setup:
            # TODO: make nodes etc.
            print(simplify, fraction)
        return bob
