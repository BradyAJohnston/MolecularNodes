import re
from pathlib import Path
from typing import TYPE_CHECKING, TypedDict, overload
import databpy as db
import numpy as np
from bpy.types import Object
from cryosparc.dataset import Dataset
from databpy.attribute import AttributeTypeNames
from scipy.spatial.transform import Rotation
from .base import Ensemble, EntityType

if TYPE_CHECKING:
    BobField = TypedDict(
        "BobField",
        {"name": str, "data": np.typing.NDArray, "atype": AttributeTypeNames},
    )

UID_RIGHT_MASK = np.uint64(2**32 - 1)
UID_REGEX = rb"(\d{21})_"


def uid_as_i32_vec(
    uids: np.ndarray[tuple[int], np.dtype[np.uint64]],
) -> np.ndarray[tuple[int, int], np.dtype[np.int32]]:
    # store the u64 UID as two i32
    left_half = (uids >> 32).astype(np.int32)
    right_half = np.bitwise_and(UID_RIGHT_MASK, uids).astype(np.int32)
    return np.column_stack((left_half, right_half))


def i32_vec_to_uid(
    split_uids: np.ndarray[tuple[int, int], np.dtype[np.int32]],
) -> np.ndarray[tuple[int], np.dtype[np.uint64]]:
    split_uids = split_uids.astype(np.uint64)
    shifts = np.zeros(split_uids.shape, dtype=np.uint32)
    shifts[:, 0] = 32
    masks = np.zeros(split_uids.shape, dtype=np.uint64)
    # have to mask the lower 32 because numpy will pad with 1s if
    # the signed int32 is negative
    masks[:, :] = 2**32 - 1
    masks = masks << shifts
    split_uids = np.bitwise_and(split_uids << shifts, masks)
    return np.bitwise_or(split_uids[:, 0], split_uids[:, 1])


class MNDataset:
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
    def shift2d(self) -> np.ndarray[tuple[int, int], np.dtype[np.float32]] | None:
        """
        Shifts aligned during 2D Classification
        """
        if (shifts := self.get("alignments2D/shift")) is not None and (
            psize := self.get("alignments2D/psize_A")
        ) is not None:
            return shifts.astype(np.float32) * psize[:, np.newaxis]

    @property
    def shift3d(self) -> np.ndarray[tuple[int, int], np.dtype[np.float32]] | None:
        """
        Aligned particle shifts in A.
        """
        if (shifts := self.get("alignments3D/shift")) is not None and (
            psize := self.get("alignments3D/psize_A")
        ) is not None:
            return shifts * psize[:, np.newaxis]

    @property
    def position(self) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
        """
        Particle positions in micrograph fractions (X, Y).
        """
        xpos = self.get("location/center_x_frac", np.zeros(len(self), dtype=np.float32))
        ypos = self.get("location/center_y_frac", np.zeros_like(xpos))
        xypos = np.column_stack((xpos, ypos))
        xypos *= (
            self.get("location/micrograph_shape", np.ones_like(xypos))
            * self.get("location/micrograph_psize_A", np.ones(len(self)))[:, np.newaxis]
        )
        zpos = np.mean(
            np.column_stack(
                [
                    self.get("ctf/df1_A", np.ones_like(xpos)),
                    self.get("ctf/df2_A", np.ones_like(xpos)),
                ]
            ),
            axis=1,
            keepdims=True,
        )
        zpos = np.median(zpos) - zpos
        return np.append(xypos, zpos, axis=1)

    @property
    def blob_path_uids(self) -> np.ndarray[tuple[int, int], np.dtype[np.int32]] | None:
        if (blob_paths := self.get("blob/path")) is None:
            return None
        try:
            # ignore type here because we handle the None.group() exception
            blob_uids = [
                np.uint64(re.search(UID_REGEX, s).group(1))  # type: ignore
                for s in blob_paths
            ]
            return uid_as_i32_vec(np.array(blob_uids))
        except (IndexError, AttributeError):
            return None

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


class CryoSPARCParticles(Ensemble):
    DEFAULT_NAME = "CryoSPARC Ensemble"

    def __init__(self, file_path: Path):
        if file_path.suffix != ".cs":
            raise ValueError(f"{file_path.name} is not a CryoSPARC .cs file")
        super().__init__(file_path=file_path)
        self._entity_type = EntityType.ENSEMBLE_CRYOSPARC
        self.dset = MNDataset(Dataset.load(self.file_path))
        if len(self.dset) == 0:
            raise ValueError("Dataset has no rows")

    def _construct_fields(self) -> "list[BobField]":
        """
        Construct named attributes from CryoSPARC fields.

        Return a sorted list of funciton argument dicts.
        We sort them here so they are presented in the same order they would be in
        a CryoSPARC dataset.
        """
        a = db.AttributeTypes
        fields: "list[BobField]" = []
        fields_to_check = {
            "ctf/df_angle_rad": a.FLOAT,
            "ctf/exp_group_id": a.INT,
            "alignments2D/shift": a.FLOAT2,
            "alignments2D/class": a.INT,
            "alignments3D/split": a.INT,
            "alignments3D/class": a.INT,
            "alignments3D/shift": a.FLOAT2,
            "alignments3D/alpha_min": a.FLOAT,
            "blob/idx": a.INT,
            "blob/shape": a.INT32_2D,
            "blob/psize_A": a.FLOAT,
            # cannot support string fields since databpy does not support strings
        }
        for name, atype in fields_to_check.items():
            if (bob_entry := self.dset.bob_entry(name=name, atype=atype)) is not None:
                fields.append(bob_entry)

        df1 = self.dset.get("ctf/df1_A")
        df2 = self.dset.get("ctf/df2_A")
        if df1 is not None and df2 is not None:
            fields.append(
                dict(
                    name="ctf/df_A",
                    data=np.mean(np.column_stack((df1, df2)), axis=1),
                    atype=a.FLOAT,
                )
            )

        loc_x = self.dset.get("location/center_x_frac")
        loc_y = self.dset.get("location/center_y_frac")
        if loc_x is not None and loc_y is not None:
            fields.append(
                dict(
                    name="location/center_frac",
                    data=np.column_stack((loc_x, loc_y)),
                    atype="FLOAT2",
                )
            )

        fields.append(
            dict(name="uid", data=uid_as_i32_vec(self.dset["uid"]), atype="INT32_2D")
        )
        if (rots := self.dset.pose3d_as_quat) is not None:
            fields.append(dict(name="alignments3D/pose", data=rots, atype="QUATERNION"))
        if (rots := self.dset.pose2d_as_quat) is not None:
            fields.append(dict(name="alignments2D/pose", data=rots, atype="QUATERNION"))
        if (blob_path_uids := self.dset.blob_path_uids) is not None:
            fields.append(
                dict(name="blob/path_uid", data=blob_path_uids, atype="INT32_2D")
            )

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
        world_scale: float = 0.1,
        fraction: float = 1.0,
        simplify=False,
    ) -> Object:
        # build positions and rotation outside _construct_fields because we want
        # them to appear first, like they would with a normal Blender object.
        self.object = db.create_object(self.dset.position * world_scale, name=name)
        self.store_named_attribute(
            data=self.dset.rotations,
            name="rotation",
            atype=db.AttributeTypes.QUATERNION,
        )
        for field in self._construct_fields():
            self.store_named_attribute(**field)
        if node_setup:
            # TODO: make nodes etc.
            print(simplify, fraction)
        return self.object
