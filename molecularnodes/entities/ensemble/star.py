from pathlib import Path
from typing import List, Optional, Union
import bpy
import databpy
import mrcfile
import numpy as np
import starfile
from databpy import AttributeTypes, BlenderObject
from pandas import CategoricalDtype, DataFrame
from PIL import Image
from scipy.spatial.transform import Rotation
from ... import blender as bl
from ...nodes import nodes
from .base import Ensemble, EntityType


class EnsembleDataFrame:
    def __init__(self, data: DataFrame) -> None:
        self.data = data
        self._coord_columns: List[str] = ["x", "y", "z"]
        self._rot_columns: List[str] = ["Rot", "Tilt", "Psi"]
        self._shift_column_names: List[str] = [
            "OriginXAngst",
            "OriginYAngst",
            "OriginZAngst",
        ]

    @property
    def coordinates(self) -> np.ndarray:
        coord = self.data[self._coord_columns].to_numpy()

        try:
            shift = self.data[self._shift_column_names].to_numpy()
            coord -= shift
        except KeyError:
            pass

        return coord

    @property
    def scale(self) -> np.ndarray:
        arr = np.zeros((len(self.data), 1), dtype=np.float32)
        arr[:] = 1.0
        return arr

    @property
    def coordinates_scaled(self) -> np.ndarray:
        return self.coordinates * self.scale

    def rotation_as_quaternion(self) -> np.ndarray:
        rot_tilt_psi_cols = self.data[self._rot_columns].to_numpy()

        quaternions = np.array(
            [
                Rotation.from_euler("ZYZ", row, degrees=True)
                .inv()
                .as_quat(scalar_first=True)
                for row in rot_tilt_psi_cols
            ]
        )
        return quaternions

    def image_id_values(self) -> np.ndarray:
        for name in [
            "rlnImageName",
            "rlnMicrographName",
            "rlnTomoName",
            "cisTEMOriginalImageFilename",
        ]:
            try:
                return self.data[name].cat.codes.to_numpy()
            except KeyError:
                pass

        return np.zeros(len(self.data), dtype=int)

    def store_data_on_object(self, obj: bpy.types.Object) -> None:
        bob = BlenderObject(obj)
        bob.store_named_attribute(
            self.rotation_as_quaternion(),
            name="rotation",
            atype=AttributeTypes.QUATERNION,
        )

        bob.store_named_attribute(
            self.image_id_values(),
            name="image_id",
            atype=AttributeTypes.INT,
        )

        for col in self.data.columns:
            if isinstance(self.data[col].dtype, CategoricalDtype):
                bob.object[f"{col}_categories"] = list(self.data[col].cat.categories)
                data = self.data[col].cat.codes.to_numpy()
                bob.store_named_attribute(data, name=col, atype=AttributeTypes.INT)
            else:
                bob.store_named_attribute(self.data[col].to_numpy(), name=col)


class RelionDataFrame(EnsembleDataFrame):
    def __init__(self, data: DataFrame) -> None:
        super().__init__(data)
        self.type: str = "relion"
        self._coord_columns = ["rlnCoordinateX", "rlnCoordinateY", "rlnCoordinateZ"]
        self._rot_columns = ["rlnAngleRot", "rlnAngleTilt", "rlnAnglePsi"]
        self._shift_column_names = [
            "rlnOriginXAngst",
            "rlnOriginYAngst",
            "rlnOriginZAngst",
        ]

    @property
    def scale(self) -> np.ndarray:
        if "rlnImagePixelSize" not in self.data:
            return super().scale

        return self.data["rlnImagePixelSize"].to_numpy().reshape((-1, 1))


class CistemDataFrame(EnsembleDataFrame):
    def __init__(self, data: DataFrame) -> None:
        super().__init__(data)
        self._adjust_defocus()
        self.type: str = "cistem"
        self._coord_columns = [
            "cisTEMOriginalXPosition",
            "cisTEMOriginalYPosition",
            "cisTEMZFromDefocus",
        ]
        self._rot_columns = ["cisTEMAnglePhi", "cisTEMAngleTheta", "cisTEMAnglePsi"]
        self._shift_column_names = [
            "origin_x",
            "origin_y",
            "origin_z",
        ]

    def _adjust_defocus(self) -> None:
        self.data["cisTEMZFromDefocus"] = (
            self.data["cisTEMDefocus1"] + self.data["cisTEMDefocus2"]
        ) / 2
        self.data["cisTEMZFromDefocus"] = (
            self.data["cisTEMZFromDefocus"] - self.data["cisTEMZFromDefocus"].median()
        )


class StarFile(Ensemble):
    object: Optional[bpy.types.Object]
    data: Optional[DataFrame]
    df: Optional[Union[RelionDataFrame, CistemDataFrame]]

    def __init__(self, file_path: Union[str, Path]) -> None:
        super().__init__(file_path)
        self.type: str = "starfile"
        self.current_image: int = -1
        self._entity_type = EntityType.ENSEMBLE_STAR

    @classmethod
    def from_starfile(cls, file_path: Union[str, Path]) -> "StarFile":
        self = cls(file_path)
        self.data = self._read()
        self.df = self._assign_df()
        return self

    @classmethod
    def from_blender_object(cls, blender_object: bpy.types.Object) -> "StarFile":
        self = cls(blender_object["starfile_path"])
        self.object = blender_object
        self.data = self._read()
        return self

    @property
    def star_node(self) -> bpy.types.Node:
        return nodes.get_star_node(self.object)

    @property
    def micrograph_material(self) -> bpy.types.Material:
        return nodes.micrograph_material()

    def _read(self) -> DataFrame:
        star_dict: dict = starfile.read(self.file_path, always_dict=True)  # type: ignore
        star: DataFrame = list(star_dict.values())[0]

        if not isinstance(star, DataFrame):
            raise ValueError("Problem opening starfile as dataframe")

        for col in star.columns:
            if star[col].dtype == "object":
                star[col] = star[col].astype("category")
        return star

    @property
    def n_images(self) -> int:
        if isinstance(self.data, dict):
            return len(self.data)
        return 1

    def _is_relion(self) -> bool:
        return (
            isinstance(self.data, dict)
            and "particles" in self.data
            and "optics" in self.data
        ) or ("rlnAnglePsi" in self.data)  # type: ignore

    def _is_cistem(self) -> bool:
        return "cisTEMAnglePsi" in self.data  # type: ignore

    def _assign_df(self) -> Union[RelionDataFrame, CistemDataFrame]:
        if self.data is None:
            raise ValueError("Data not loaded. Call from_starfile() first.")
        if self._is_relion():
            return RelionDataFrame(self.data)
        elif self._is_cistem():
            return CistemDataFrame(self.data)
        else:
            raise ValueError(
                "File is not a valid RELION>=3.1 or cisTEM STAR file, other formats are not currently supported."
            )

    def _convert_mrc_to_tiff(self) -> Path:
        if self.object is None:
            raise ValueError("Object not set. Call from_blender_object() first.")

        if self._is_relion():
            micrograph_path = self.object["rlnMicrographName_categories"][
                self.star_node.inputs["Image"].default_value - 1  # type: ignore
            ]
        elif self._is_cistem():
            micrograph_path = self.object["cisTEMOriginalImageFilename_categories"][
                self.star_node.inputs["Image"].default_value - 1  # type: ignore
            ].strip("'")
        else:
            raise ValueError("File is not a valid RELION>=3.1 or cisTEM STAR file")

        micrograph_path = Path(micrograph_path)
        if not micrograph_path.exists():
            pot_micrograph_path = Path(self.file_path).parent / micrograph_path
            if not pot_micrograph_path.exists():
                if self._is_relion():
                    pot_micrograph_path = (
                        Path(self.file_path).parent.parent.parent / micrograph_path
                    )
                    if not pot_micrograph_path.exists():
                        raise FileNotFoundError(
                            f"Micrograph file {micrograph_path} not found"
                        )
                else:
                    raise FileNotFoundError(
                        f"Micrograph file {micrograph_path} not found"
                    )
            micrograph_path = pot_micrograph_path

        tiff_path = Path(micrograph_path).with_suffix(".tiff")
        if not tiff_path.exists():
            with mrcfile.open(micrograph_path) as mrc:
                if mrc.data is None:
                    raise ValueError(f"Micrograph file {micrograph_path} is empty")
                micrograph_data = mrc.data.copy()

            if micrograph_data.ndim == 3:
                micrograph_data = np.sum(micrograph_data, axis=0)
            micrograph_data = (micrograph_data - micrograph_data.min()) / (
                micrograph_data.max() - micrograph_data.min()
            )

            if micrograph_data.dtype != np.float32:
                micrograph_data = micrograph_data.astype(np.float32)

            Image.fromarray(micrograph_data[::-1, :]).save(tiff_path)
        return tiff_path

    def create_object(
        self,
        name: str = "StarFileObject",
        node_setup: bool = True,
        world_scale: float = 0.01,
        fraction: float = 1.0,
        simplify: bool = True,
    ) -> bpy.types.Object:
        if self.df is None:
            raise ValueError("DataFrame not assigned. Call from_starfile() first.")

        self.object = databpy.create_object(
            self.df.coordinates_scaled * world_scale, collection=bl.coll.mn(), name=name
        )
        self.object.mn.entity_type = self._entity_type.value
        self.df.store_data_on_object(self.object)

        if node_setup:
            nodes.create_starting_nodes_starfile(self.object)

        self.object["starfile_path"] = str(self.file_path)
        return self.object
