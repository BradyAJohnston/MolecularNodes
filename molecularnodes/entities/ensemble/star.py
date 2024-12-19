from pathlib import Path

import bpy
import mrcfile
import numpy as np
import starfile
from PIL import Image

from ... import blender as bl
from databpy import AttributeTypes
import databpy
from .ensemble import Ensemble


class StarFile(Ensemble):
    def __init__(self, file_path):
        super().__init__(file_path)
        self.type = "starfile"

    @classmethod
    def from_starfile(cls, file_path):
        self = cls(file_path)
        self.data = self._read()
        self.star_type = None
        self.positions = None
        self.current_image = -1
        self._create_mn_columns()
        self.n_images = self._n_images()
        return self

    @classmethod
    def from_blender_object(cls, blender_object):
        self = cls(blender_object["starfile_path"])
        self.object = blender_object
        self.data = self._read()
        self.star_type = None
        self.positions = None
        self.current_image = -1
        self._create_mn_columns()
        self.n_images = self._n_images()
        bpy.app.handlers.depsgraph_update_post.append(self._update_micrograph_texture)
        return self

    @property
    def star_node(self):
        return bl.nodes.get_star_node(self.object)

    @property
    def micrograph_material(self):
        return bl.nodes.MN_micrograph_material()

    def _read(self):
        star = starfile.read(self.file_path)
        return star

    def _n_images(self):
        if isinstance(self.data, dict):
            return len(self.data)
        return 1

    def _create_mn_columns(self):
        # only RELION 3.1 and cisTEM STAR files are currently supported, fail gracefully
        if (
            isinstance(self.data, dict)
            and "particles" in self.data
            and "optics" in self.data
        ):
            self.star_type = "relion"
        elif "cisTEMAnglePsi" in self.data:
            self.star_type = "cistem"
        else:
            raise ValueError(
                "File is not a valid RELION>=3.1 or cisTEM STAR file, other formats are not currently supported."
            )

        # Get absolute position and orientations
        if self.star_type == "relion":
            df = self.data["particles"].merge(self.data["optics"], on="rlnOpticsGroup")

            # get necessary info from dataframes
            # Standard cryoEM starfile don't have rlnCoordinateZ. If this column is not present
            # Set it to "0"
            if "rlnCoordinateZ" not in df:
                df["rlnCoordinateZ"] = 0

            self.positions = df[
                ["rlnCoordinateX", "rlnCoordinateY", "rlnCoordinateZ"]
            ].to_numpy()
            pixel_size = df["rlnImagePixelSize"].to_numpy().reshape((-1, 1))
            self.positions = self.positions * pixel_size
            shift_column_names = [
                "rlnOriginXAngst",
                "rlnOriginYAngst",
                "rlnOriginZAngst",
            ]
            if all([col in df.columns for col in shift_column_names]):
                shifts_ang = df[shift_column_names].to_numpy()
                self.positions -= shifts_ang
            df["MNAnglePhi"] = df["rlnAngleRot"]
            df["MNAngleTheta"] = df["rlnAngleTilt"]
            df["MNAnglePsi"] = df["rlnAnglePsi"]
            df["MNPixelSize"] = df["rlnImagePixelSize"]
            try:
                df["MNImageId"] = (
                    df["rlnMicrographName"].astype("category").cat.codes.to_numpy()
                )
            except KeyError:
                try:
                    df["MNImageId"] = (
                        df["rlnTomoName"].astype("category").cat.codes.to_numpy()
                    )
                except KeyError:
                    df["MNImageId"] = 0.0

            self.data = df

        elif self.star_type == "cistem":
            df = self.data
            df["cisTEMZFromDefocus"] = (df["cisTEMDefocus1"] + df["cisTEMDefocus2"]) / 2
            df["cisTEMZFromDefocus"] = (
                df["cisTEMZFromDefocus"] - df["cisTEMZFromDefocus"].median()
            )
            self.positions = df[
                [
                    "cisTEMOriginalXPosition",
                    "cisTEMOriginalYPosition",
                    "cisTEMZFromDefocus",
                ]
            ].to_numpy()
            df["MNAnglePhi"] = df["cisTEMAnglePhi"]
            df["MNAngleTheta"] = df["cisTEMAngleTheta"]
            df["MNAnglePsi"] = df["cisTEMAnglePsi"]
            df["MNPixelSize"] = df["cisTEMPixelSize"]
            df["MNImageId"] = (
                df["cisTEMOriginalImageFilename"]
                .astype("category")
                .cat.codes.to_numpy()
            )

    def _convert_mrc_to_tiff(self):
        if self.star_type == "relion":
            micrograph_path = self.object["rlnMicrographName_categories"][
                self.star_node.inputs["Image"].default_value - 1
            ]
        elif self.star_type == "cistem":
            micrograph_path = self.object["cisTEMOriginalImageFilename_categories"][
                self.star_node.inputs["Image"].default_value - 1
            ].strip("'")
        else:
            return False

        # This could be more elegant
        if not Path(micrograph_path).exists():
            pot_micrograph_path = Path(self.file_path).parent / micrograph_path
            if not pot_micrograph_path.exists():
                if self.star_type == "relion":
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
                micrograph_data = mrc.data.copy()

            # For 3D data sum over the z axis. Probalby would be nicer to load the data as a volume
            if micrograph_data.ndim == 3:
                micrograph_data = np.sum(micrograph_data, axis=0)
            # Normalize the data to 0-1
            micrograph_data = (micrograph_data - micrograph_data.min()) / (
                micrograph_data.max() - micrograph_data.min()
            )

            if micrograph_data.dtype != np.float32:
                micrograph_data = micrograph_data.astype(np.float32)

            # Need to invert in Y to generate the correct tiff
            Image.fromarray(micrograph_data[::-1, :]).save(tiff_path)
        return tiff_path

    def _update_micrograph_texture(self, *_):
        try:
            show_micrograph = self.star_node.inputs["Show Micrograph"]
            _ = self.object["mn"]
        except ReferenceError:
            bpy.app.handlers.depsgraph_update_post.remove(
                self._update_micrograph_texture
            )
            return
        if self.star_node.inputs["Image"].default_value == self.current_image:
            return
        else:
            self.current_image = self.star_node.inputs["Image"].default_value
        if not show_micrograph:
            return
        tiff_path = self._convert_mrc_to_tiff()
        if tiff_path:
            try:
                image_obj = bpy.data.images[tiff_path.name]
            except KeyError:
                image_obj = bpy.data.images.load(str(tiff_path))
            image_obj.colorspace_settings.name = "Non-Color"
            self.micrograph_material.node_tree.nodes["Image Texture"].image = image_obj
            self.star_node.inputs["Micrograph"].default_value = image_obj

    def create_object(self, name="StarFileObject", node_setup=True, world_scale=0.01):
        self.object = databpy.create_object(
            self.positions * world_scale, collection=bl.coll.mn(), name=name
        )

        self.object.mn["entity_type"] = "star"

        # create attribute for every column in the STAR file
        for col in self.data.columns:
            col_type = self.data[col].dtype
            # If col_type is numeric directly add
            if np.issubdtype(col_type, np.number):
                self.store_named_attribute(
                    name=col,
                    data=self.data[col].to_numpy().reshape(-1),
                    atype=AttributeTypes.FLOAT,
                )

            # If col_type is object, convert to category and add integer values
            elif isinstance(col_type, object):
                codes = (
                    self.data[col].astype("category").cat.codes.to_numpy().reshape(-1)
                )
                self.store_named_attribute(
                    data=codes, name=col, atype=AttributeTypes.INT
                )
                # Add the category names as a property to the blender object
                self.object[f"{col}_categories"] = list(
                    self.data[col].astype("category").cat.categories
                )

        if node_setup:
            bl.nodes.create_starting_nodes_starfile(self.object, n_images=self.n_images)

        self.object["starfile_path"] = str(self.file_path)
        bpy.app.handlers.depsgraph_update_post.append(self._update_micrograph_texture)
        return self.object
