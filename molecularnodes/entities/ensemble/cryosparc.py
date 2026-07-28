from pathlib import Path
from typing import TYPE_CHECKING, TypedDict, overload
import databpy as db
import numpy as np
from databpy.attribute import AttributeTypeNames
from scipy.spatial.transform import Rotation
from .base import Ensemble, EntityType

if TYPE_CHECKING:
    from molecularnodes.ui.ops import MN_OT_Import_CryoSPARC_File

    BobField = TypedDict(
        "BobField",
        {"name": str, "data": np.typing.NDArray, "atype": AttributeTypeNames},
    )


class Dataset:
    def __init__(self, dset: np.typing.NDArray, op: "MN_OT_Import_CryoSPARC_File"):
        self.dset = dset
        self.op = op

    def __getitem__(self, key: str) -> np.typing.NDArray:
        return self.dset[key]

    def __contains__(self, key: str) -> bool:
        if self.dset.dtype.names is None:
            return False
        return key in self.dset.dtype.names

    def __len__(self) -> int:
        return len(self.dset)

    def bob_entry(self, name, atype) -> "BobField | None":
        if (data := self.get(name)) is None:
            return None
        else:
            return dict(data=data, name=name, atype=atype)

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
    def uid_as_i32(self) -> np.ndarray[tuple[int], np.dtype[np.int32]]:
        uids = self["uid"].astype(np.int32)
        counts = np.unique(uids, return_counts=True)[1]
        if np.any(counts != 1):
            self.op.report({"WARNING"}, "UID collisions due to i32 coercion.")
        return uids

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
    def pose2d_as_quat(
        self,
    ) -> np.ndarray[tuple[int, int], np.dtype[np.float32]] | None:
        if (poses := self.get("alignments2D/pose")) is None:
            return None
        converted = np.zeros((len(self), 3), dtype=np.float32)
        converted[:, 2] = poses
        return Rotation.from_rotvec(converted).as_quat(scalar_first=True)

    @property
    def pose3d_as_quat(
        self,
    ) -> np.ndarray[tuple[int, int], np.dtype[np.float32]] | None:
        if (poses := self.get("alignments3D/pose")) is None:
            return None
        return Rotation.from_rotvec(poses).as_quat(scalar_first=True)

    @property
    def rotations(self) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
        """
        Pose 3D
        """
        if (rots := self.pose3d_as_quat) is not None:
            return rots
        if (rots := self.pose2d_as_quat) is not None:
            return rots
        return np.zeros((len(self), 4), dtype=np.float32)


class CryoSPARC(Ensemble):
    DEFAULT_NAME = "CryoSPARC Ensemble"

    def __init__(self, file_path: Path, op: "MN_OT_Import_CryoSPARC_File"):
        if file_path.suffix != ".cs":
            raise ValueError(f"{file_path.name} is not a CryoSPARC .cs file")
        super().__init__(file_path=file_path)
        self._entity_type = EntityType.ENSEMBLE_CRYOSPARC
        self.op = op
        self.dset = Dataset(np.load(self.file_path), op=op)
        if len(self.dset) == 0:
            raise ValueError("Dataset has no rows")

    def _construct_fields(self) -> "list[BobField]":
        """
        Construct named attributes from CryoSPARC fields.

        Return a sorted list of funciton argument dicts.
        We sort them here so they are presented in the same order they would be in
        a CryoSPARC dataset.
        """
        fields: "list[BobField]" = []
        if not np.allclose((defocus := self.dset.defocus), 0.0):
            fields.append(dict(data=defocus[:, 0], name="ctf/df1_um", atype="FLOAT"))
            fields.append(dict(data=defocus[:, 1], name="ctf/df2_um", atype="FLOAT"))
        fields_to_check = {
            "ctf/df_angle_rad": "FLOAT",
            "ctf/exp_group_id": "INT",
            "alignments2D/shift": "FLOAT2",
            "alignments2D/class": "INT",
            "alignments3D/split": "INT",
            "alignments3D/class": "INT",
            "alignments3D/shift": "FLOAT2",
            "alignments3D/alpha_min": "FLOAT",
            # cannot support string fields since databpy does not support strings
        }
        for name, atype in fields_to_check.items():
            if (bob_entry := self.dset.bob_entry(name=name, atype=atype)) is not None:
                fields.append(bob_entry)

        fields.append(dict(name="uid", data=self.dset.uid_as_i32, atype="INT"))
        if (rots := self.dset.pose3d_as_quat) is not None:
            fields.append(dict(name="alignments3D/pose", data=rots, atype="QUATERNION"))
        if (rots := self.dset.pose2d_as_quat) is not None:
            fields.append(dict(name="alignments2D/pose", data=rots, atype="QUATERNION"))

        component_coord = 0
        while (
            component_val := self.dset.get(f"components_mode_{component_coord}/value")
        ) is not None:
            fields.append(
                dict(
                    name=f"components_mode_{component_coord}/value",
                    data=component_val,
                    atype="FLOAT",
                )
            )
            component_coord += 1

        return sorted(fields, key=lambda field: field.get("name"))

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
        position_scaler = np.ones_like(self.dset.positions)
        position_scaler[:, 2] = world_scale
        bob = db.create_bob(self.dset.positions * position_scaler, name=name)
        bob.store_named_attribute(
            data=self.dset.rotations, name="rotation", atype="QUATERNION"
        )
        for field in self._construct_fields():
            bob.store_named_attribute(**field)
        if node_setup:
            # TODO: make nodes etc.
            print(simplify, fraction)
        return bob
