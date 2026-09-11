import json
from pathlib import Path
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
from ...nodes import geometry
from ...nodes.utils import get_star_node
from .base import Ensemble, EntityType


class EnsembleDataFrame:
    def __init__(self, data: DataFrame) -> None:
        self._world_scale = 0.1
        self.data = data
        self._coord_columns: list[str] = ["x", "y", "z"]
        self._rot_columns: list[str] = ["Rot", "Tilt", "Psi"]
        self._shift_column_names: list[str] = [
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
        return self.coordinates * self.scale * self._world_scale

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
            self.image_id_values(),
            name="image_id",
            atype=AttributeTypes.INT,
        )

        categories: dict[str, list] = {}
        for col in self.data.columns:
            if isinstance(self.data[col].dtype, CategoricalDtype):
                categories[col] = self.data[col].cat.categories.tolist()
                data = self.data[col].cat.codes.to_numpy()
                bob.store_named_attribute(data, name=col, atype=AttributeTypes.INT)
            else:
                bob.store_named_attribute(self.data[col].to_numpy(), name=col)

        bob.object.mn.categories = categories


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


class NDJSONDataFrame(EnsembleDataFrame):
    """
    Point annotations from a CZI CryoET Data Portal ``.ndjson`` file.

    Each line is a JSON object with a ``location`` and, for oriented points, a
    3x3 ``xyz_rotation_matrix`` that rotates the reference structure into the
    tomogram frame. Coordinates are in voxels of the parent tomogram - the file
    carries no pixel size, so no additional scaling is applied on import.
    """

    _matrix_columns = [f"rotation_{i}{j}" for i in range(3) for j in range(3)]

    def __init__(self, data: DataFrame) -> None:
        super().__init__(data)
        self.type: str = "ndjson"

    @property
    def _has_rotation(self) -> bool:
        return all(column in self.data for column in self._matrix_columns)

    def transforms(self) -> np.ndarray:
        """
        The full 4x4 world-scaled transform for each point.

        The rotation part is the file's ``xyz_rotation_matrix`` (identity for
        plain points) and the translation is the world-scaled position.
        """
        n_points = len(self.data)
        transforms = np.tile(np.identity(4, dtype=float), (n_points, 1, 1))
        if self._has_rotation:
            transforms[:, :3, :3] = (
                self.data[self._matrix_columns].to_numpy().reshape(n_points, 3, 3)
            )
        transforms[:, :3, 3] = self.coordinates_scaled
        return transforms

    def store_data_on_object(self, obj: bpy.types.Object) -> None:
        bob = BlenderObject(obj)

        bob.store_named_attribute(
            self.image_id_values(),
            name="image_id",
            atype=AttributeTypes.INT,
        )
        # Blender stores float4x4 attributes column-major, so the row-major
        # numpy matrices must be transposed or GN sees the inverse rotation
        bob.store_named_attribute(
            self.transforms().transpose(0, 2, 1),
            name="transform",
            atype=AttributeTypes.FLOAT4X4,
        )
        if "instance_id" in self.data:
            bob.store_named_attribute(
                self.data["instance_id"].to_numpy(dtype=int),
                name="instance_id",
                atype=AttributeTypes.INT,
            )


class StarFile(Ensemble):
    data_reader: DataFrame | None
    data_frame: RelionDataFrame | CistemDataFrame | NDJSONDataFrame | None

    def __init__(self, file_path: str | Path) -> None:
        super().__init__(file_path)
        self.type: str = "starfile"
        self.current_image: int = -1
        self._entity_type = EntityType.ENSEMBLE_STAR
        self.data_reader = self._read()
        self.data_frame = self._assign_df()

    @classmethod
    def from_blender_object(cls, blender_object: bpy.types.Object) -> "StarFile":
        self = cls(blender_object.mn.filepath)
        self.object = blender_object
        return self

    @property
    def star_node(self) -> bpy.types.Node:
        return get_star_node(self.object)

    @property
    def micrograph_material(self) -> bpy.types.Material:
        return bpy.data.materials["MN_micrograph_material"]

    def _read(self) -> DataFrame:
        if self._is_ndjson():
            return self._read_ndjson()

        star_dict: dict = starfile.read(self.file_path, always_dict=True)  # type: ignore
        star: DataFrame = list(star_dict.values())[0]

        if not isinstance(star, DataFrame):
            raise ValueError("Problem opening starfile as dataframe")

        for col in star.columns:
            if star[col].dtype == "object":
                star[col] = star[col].astype("category")
        return star

    def _read_ndjson(self) -> DataFrame:
        """
        Read a CZI CryoET Data Portal ``.ndjson`` annotation file.

        Handles ``point``, ``orientedPoint`` and ``instancePoint`` annotations:
        one JSON object per line with a ``location``, and optionally a 3x3
        ``xyz_rotation_matrix`` (stored element-wise in ``rotation_{ij}``
        columns) and an ``instance_id``.
        """
        with open(self.file_path) as file:
            records: list[dict] = [json.loads(line) for line in file if line.strip()]

        data = DataFrame(
            {axis: [record["location"][axis] for record in records] for axis in "xyz"}
        )

        matrices = [record.get("xyz_rotation_matrix") for record in records]
        if any(matrix is not None for matrix in matrices):
            stacked = np.array(
                [
                    matrix if matrix is not None else np.identity(3)
                    for matrix in matrices
                ],
                dtype=float,
            )
            for i in range(3):
                for j in range(3):
                    data[f"rotation_{i}{j}"] = stacked[:, i, j]

        instance_ids = [record.get("instance_id") for record in records]
        if all(instance_id is not None for instance_id in instance_ids):
            data["instance_id"] = np.array(instance_ids, dtype=int)

        return data

    @property
    def n_images(self) -> int:
        if isinstance(self.data_reader, dict):
            return len(self.data_reader)
        return 1

    def _is_ndjson(self) -> bool:
        return self.file_path.suffix == ".ndjson"

    def _is_relion(self) -> bool:
        return (
            isinstance(self.data_reader, dict)
            and "particles" in self.data_reader
            and "optics" in self.data_reader
        ) or ("rlnAnglePsi" in self.data_reader)  # type: ignore

    def _is_cistem(self) -> bool:
        return "cisTEMAnglePsi" in self.data_reader  # type: ignore

    def _assign_df(self) -> RelionDataFrame | CistemDataFrame | NDJSONDataFrame:
        if self.data_reader is None:
            raise ValueError("Data not loaded. Call StarFile.load() first.")
        if self._is_ndjson():
            return NDJSONDataFrame(self.data_reader)
        elif self._is_relion():
            return RelionDataFrame(self.data_reader)
        elif self._is_cistem():
            return CistemDataFrame(self.data_reader)
        else:
            raise ValueError(
                "File is not a valid RELION>=3.1 or cisTEM STAR file, other formats are not currently supported."
            )

    def _convert_mrc_to_tiff(self) -> Path:
        if self.object is None:
            raise ValueError("Object not set. Call from_blender_object() first.")

        categories = self.props.categories
        image_index = self.star_node.inputs["Image"].default_value - 1
        if self._is_relion():
            micrograph_path = categories["rlnMicrographName"][image_index]
        elif self._is_cistem():
            micrograph_path = categories["cisTEMOriginalImageFilename"][
                image_index
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
    ) -> bpy.types.Object:
        if self.data_frame is None:
            raise ValueError("DataFrame not assigned. Call StarFile.load() first.")

        self.object = databpy.create_object(
            self.data_frame.coordinates_scaled,
            collection=bl.coll.mn(),
            name=name,
        )
        self.props.entity_type = self._entity_type.value
        self.data_frame.store_data_on_object(self.object)

        if node_setup:
            with self.tree.reset() as (input, join):
                input >> geometry.StarfileInstances() >> join

        self.props.filepath = str(self.file_path)
        return self.object
