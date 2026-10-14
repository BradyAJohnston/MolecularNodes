"""
Import a `bentopy <https://github.com/marrink-lab/bentopy>`_ placement list.

``bentopy-pack`` writes a JSON *placement list* that records, for each packed
segment, the structure file it came from and batches of positions that share a
rotation. ``bentopy-render`` later expands this into a (potentially huge) ``.gro``
file. Here the placement list is imported directly as an ensemble: each segment's
structure is loaded once and instanced at every placement, so large packed systems
stay light in Blender.
"""

import json
from pathlib import Path
from typing import cast
import bpy
import numpy as np
from biotite.structure import AtomArray, AtomArrayStack
from biotite.structure.io import gro
from nodebpy.builder import TreeBuilder
from ... import blender as bl
from ... import material
from ...nodes import geometry
from ..molecule.pdb import PDBReader
from ..molecule.reader import ReaderBase
from ..utilities import create_object
from .base import Ensemble, EntityType


def _read_bentopy_coords(path: Path) -> np.ndarray:
    """
    Read atom positions (nm) exactly as ``bentopy-render`` does.

    The placement transform depends on the mean position and minimum corner of
    the structure as bentopy parsed it, so this mirrors its parser (every
    ``ATOM``/``HETATM`` record of a PDB, every atom line of a GRO) rather than
    the biotite reader used for display, which may drop alternate locations.
    """
    lines = path.read_text().splitlines()
    if path.suffix == ".gro":
        n_atoms = int(lines[1].strip())
        atoms = lines[2 : 2 + n_atoms]
        return np.array(
            [[float(a[20:28]), float(a[28:36]), float(a[36:44])] for a in atoms]
        )
    atoms = [a for a in lines if a.startswith(("ATOM", "HETATM"))]
    coords = [[float(a[30:38]), float(a[38:46]), float(a[46:54])] for a in atoms]
    # PDB coordinates are in Å
    return np.array(coords) / 10


def _read_structure(path: Path) -> AtomArray:
    if path.suffix == ".gro":
        array = gro.GROFile.read(path).get_structure(model=1)
        array = ReaderBase.set_standard_annotations(array)
    else:
        array = PDBReader(path).array
    if isinstance(array, AtomArrayStack):
        array = array[0]
    return array


class BentopyReader:
    """
    Parse a bentopy placement list into per-instance transforms.

    Parameters
    ----------
    file_path : str | Path
        Path to the placement list ``.json`` written by ``bentopy-pack``.
    root : str | Path | None, optional
        Directory that relative structure paths in the placement list are resolved
        against, like ``bentopy-render --root``. Defaults to the directory
        containing the placement list.
    """

    def __init__(self, file_path: str | Path, root: str | Path | None = None):
        self.file_path = Path(file_path)
        self.root = Path(root) if root is not None else self.file_path.parent
        with open(self.file_path) as f:
            self.data: dict = json.load(f)
        self.title: str = self.data.get("title", "")
        self.size = np.array(self.data["size"], dtype=float)
        self.placements: list[dict] = self.data["placements"]

    def structure_path(self, placement: dict) -> Path:
        path = Path(placement["path"])
        return path if path.is_absolute() else self.root / path

    @property
    def names(self) -> list[str]:
        return [p["name"] for p in self.placements]

    def get_transforms(self) -> np.ndarray:
        """
        The transform of every placed instance, in placement-list order.

        ``bentopy-render`` centres each structure on its mean position ``c``,
        rotates it by the batch rotation ``R`` and then places the minimum corner
        ``m`` of the rotated structure's bounding box (clamped to include the
        origin) at the listed position ``p``. For an untransformed structure ``x``
        that is ``R @ (x - c) - m + p``, so the instance transform is ``R`` with a
        translation of ``p - m - R @ c``.

        Returns
        -------
        np.ndarray
            Structured array with a ``(4, 4)`` ``transform`` (translation in Å)
            and the ``chain_id`` (placement index) of each instance.
        """
        dtype = [("chain_id", int), ("transform", float, (4, 4))]
        transforms = []
        for i, placement in enumerate(self.placements):
            coords = _read_bentopy_coords(self.structure_path(placement))
            centre = coords.mean(axis=0)
            for batch in placement["batches"]:
                positions = np.array(batch["positions"], dtype=float).reshape(-1, 3)
                if len(positions) == 0:
                    continue
                # rotations are stored row-major
                rotation = np.array(batch["rotation"], dtype=float)
                rotated = (coords - centre) @ rotation.T
                corner = np.minimum(rotated.min(axis=0), 0.0)
                arr = np.zeros(len(positions), dtype=dtype)
                arr["chain_id"] = i
                arr["transform"] = np.eye(4)
                arr["transform"][:, :3, :3] = rotation
                # nm -> Å, which `create_data_object` then scales to Blender units
                translation = positions - corner - rotation @ centre
                arr["transform"][:, :3, 3] = translation * 10
                transforms.append(arr)

        if not transforms:
            return np.zeros(0, dtype=dtype)
        return np.concatenate(transforms)


class Bentopy(Ensemble):
    """
    A packed system from a bentopy placement list, imported as instances.

    Parameters
    ----------
    file_path : str | Path
        Path to the placement list ``.json`` written by ``bentopy-pack``.
    root : str | Path | None, optional
        Directory that relative structure paths are resolved against. Defaults to
        the directory containing the placement list.
    """

    def __init__(self, file_path: str | Path, root: str | Path | None = None) -> None:
        super().__init__(file_path)
        self._entity_type = EntityType.ENSEMBLE_BENTOPY
        self.file = BentopyReader(self.file_path, root=root)
        self.transformations = self.file.get_transforms()

    def create_object(
        self,
        name: str = "Bentopy",
        node_setup: bool = True,
    ) -> bpy.types.Object:
        self._create_data_object(name=name)
        # Blender may have renamed the object to avoid a clash (`name.001`), so
        # name the instances after the object rather than the requested name
        collection = self._create_object_instances(
            name=self.object.name, node_setup=node_setup
        )

        with self.tree.reset() as (atoms, join):
            (atoms >> geometry.EnsembleInstance(instances=collection) >> join)

        return self.object

    def _create_object_instances(
        self, name: str = "Bentopy", node_setup: bool = True
    ) -> bpy.types.Collection:
        collection = bl.coll.bentopy(name)

        if node_setup:
            # a new group rather than `TreeBuilder.geometry(name)`, which would
            # build into an existing group of the same name from an earlier import
            tree = bpy.data.node_groups.new(
                f"MN_pack_instance_{name}", "GeometryNodeTree"
            )
            with TreeBuilder.geometry(tree) as instance_tree:
                atoms = instance_tree.inputs.geometry("Atoms")
                style = geometry.StyleSpheres(
                    material=material.AmbientOcclusion().material
                )

                atoms >> style >> instance_tree.outputs.geometry("Geometry")
            instance_tree.tree.is_modifier = True

        for i, placement in enumerate(self.file.placements):
            array = _read_structure(self.file.structure_path(placement))
            # Collection Info sorts separated children alphabetically, so the index
            # prefix keeps each object at the position its `chain_id` points to
            obj = create_object(
                array=array,
                name=f"{i:04d}_{placement['name']}",
                collection=collection,
            )
            obj.mn.entity_type = self._entity_type.value

            if node_setup:
                mod = cast(
                    bpy.types.NodesModifier,
                    obj.modifiers.new("Molecular Nodes", "NODES"),
                )
                mod.node_group = instance_tree.tree

        self.instance_collection = collection

        return collection

    def _create_data_object(self, name: str = "DataObject") -> None:
        self.object = bl.mesh.create_data_object(
            self.transformations, name=name, collection=bl.coll.mn()
        )
        self.props.chain_ids = self.file.names
        self.props.entity_type = self._entity_type.value
