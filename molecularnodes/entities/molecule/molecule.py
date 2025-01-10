import io
import time
import warnings
from abc import ABCMeta
from pathlib import Path
from typing import Optional, Tuple, Union
import json

import biotite.structure as struc
import bpy
import numpy as np
import numpy.typing as npt
from biotite import InvalidFileError

from ... import blender as bl
from ... import color, data, utils
from databpy import Domains, AttributeTypes
import databpy
from ..entity import MolecularEntity, EntityType


class Molecule(MolecularEntity, metaclass=ABCMeta):
    """
    Abstract base class for representing a molecule.

    It associates the atomic data (the array) with the created 3D model inside of Blender
    (the object). If multiple conformations are imported, then a `frames` collection
    is also instantiated.

    The `named_attribute()` and `store_named_attribute()` methods access and set attributes on
    `object` that is in the Blender scene.

    Attributes
    ----------
    file_path : str
        The file path to the file which stores the atomic coordinates.
    file : Any
        The opened file.
    object : bpy.types.Object
        The Blender object representing the molecule.
    frames : bpy.types.Collection
        The Blender collection which holds the objects making up the frames to animate.
    array: np.ndarray:
        The numpy array which stores the atomic coordinates and associated attributes.
    entity_ids : np.ndarray
        The entity IDs of the molecule.
    chain_ids : np.ndarray
        The chain IDs of the molecule.

    Methods
    -------
    store_named_attribute(data, name='NewAttribute', type=None, domain='POINT', overwrite=True)
        Set an attribute on the object for the molecule.
    named_attribute(name='position')
        Get the value of an attribute on the object for the molecule.
    create_object(name='NewMolecule', style='spheres', selection=None, build_assembly=False, centre='', del_solvent=True, collection=None, verbose=False)
        Create a 3D model for the molecule, based on the values from self.array.
    assemblies(as_array=False)
        Get the biological assemblies of the molecule.
    """

    def __init__(self, file_path: Union[str, Path, io.BytesIO]):
        """
        Initialize the Molecule object.

        Parameters
        ----------
        file_path : Union[str, Path, io.BytesIO]
            The file path to the file which stores the atomic coordinates.
        """
        super().__init__()
        self._parse_filepath(file_path=file_path)
        self.file: str
        self.array: np.ndarray
        self._frames_collection: str | None
        self._entity_type = EntityType.MOLECULE

    @property
    def frames(self) -> bpy.types.Collection:
        """
        Get the collection of frames for the molecule.

        Returns
        -------
        bpy.types.Collection
            The collection of frames for the molecule.
        """
        if self._frames_collection is None:
            raise ValueError("No frames collection has been set for this molecule.")

        return bpy.data.collections[self._frames_collection]

    @frames.setter
    def frames(self, value: bpy.types.Collection):
        """
        Set the collection of frames for the molecule.

        Parameters
        ----------
        value : bpy.types.Collection
            The collection of frames for the molecule.
        """
        if value is None:
            self._frames_collection = None
            return
        if not isinstance(value, bpy.types.Collection):
            raise TypeError("The frames must be a bpy.types.Collection.")

        self._frames_collection = value.name

    @classmethod
    def _read(self, file_path: Union[Path, io.BytesIO]):
        """
        Initially open the file, ready to extract the required data.

        Parameters
        ----------
        file_path : Union[Path, io.BytesIO]
            The file path to the file which stores the atomic coordinates.
        """
        pass

    def _parse_filepath(self, file_path: Union[Path, str, io.BytesIO]) -> None:
        """
        If this is an actual file resolve the path - if a bytes IO resolve this as well.

        Parameters
        ----------
        file_path : Union[Path, str, io.BytesIO]
            The file path to the file which stores the atomic coordinates.
        """
        if isinstance(file_path, io.BytesIO):
            self.file = self._read(file_path=file_path)
        elif isinstance(file_path, io.StringIO):
            self.file = self._read(file_path=file_path)
        else:
            self.file_path = bl.path_resolve(file_path)
            self.file = self._read(self.file_path)

    def __len__(self) -> int:
        """
        Get the number of atoms in the molecule.

        Returns
        -------
        int
            The number of atoms in the molecule.
        """
        if hasattr(self, "object"):
            if self.object:
                return len(self.object.data.vertices)
        if self.array:
            return len(self.array)
        else:
            return 0

    @property
    def n_models(self):
        """
        Get the number of models in the molecule.

        Returns
        -------
        int
            The number of models in the molecule.
        """
        if isinstance(self.array, struc.AtomArray):
            return 1
        else:
            return self.array.shape[0]

    @property
    def tree(self) -> bpy.types.GeometryNodeTree:
        return self.object.modifiers["MolecularNodes"].node_group

    @property
    def chain_ids(self) -> Optional[list]:
        """
        Get the unique chain IDs of the molecule.

        Returns
        -------
        Optional[list]
            The unique chain IDs of the molecule, or None if not available.
        """
        if self.array:
            if hasattr(self.array, "chain_id"):
                return np.unique(self.array.chain_id).tolist()

        return None

    def create_object(
        self,
        name: str = "NewMolecule",
        style: str = "spheres",
        selection: np.ndarray = None,
        build_assembly=False,
        centre: str = "",
        del_solvent: bool = True,
        del_hydrogen: bool = False,
        collection=None,
        verbose: bool = False,
        color: Optional[str] = "common",
    ) -> bpy.types.Object:
        """
        Create a 3D model of the molecule inside of Blender.

        Creates a 3D model with one vertex per atom, and one edge per bond. Each vertex
        is given attributes which correspond to the atomic data such as `atomic_number` for
        the element and `res_name` for the residue name that the atom is associated with.

        If multiple conformations of the structure are detected, the collection attribute
        is also created which will store an object for each conformation, so that the
        object can interpolate between those conformations.

        Parameters
        ----------
        name : str, optional
            The name of the model. Default is 'NewMolecule'.
        style : str, optional
            The style of the model. Default is 'spheres'.
        selection : np.ndarray, optional
            The selection of atoms to include in the model. Default is None.
        build_assembly : bool, optional
            Whether to build the biological assembly. Default is False.
        centre : str, optional
            Denote method used to determine center of structure. Default is '',
            resulting in no translational motion being removed. Accepted values
            are `centroid` or `mass`. Any other value will result in default
            behavior.
        del_solvent : bool, optional
            Whether to delete solvent molecules. Default is True.
        del_hydrogen: bool, optional
            Whether to delete hydrogen atoms. Default is False.
        collection : str, optional
            The collection to add the model to. Default is None.
        verbose : bool, optional
            Whether to print verbose output. Default is False.
        color : Optional[str], optional
            The color scheme to use for the model. Default is 'common'.

        Returns
        -------
        bpy.types.Object
            The created 3D model, as an object in the 3D scene.
        """
        is_stack = isinstance(self.array, struc.AtomArrayStack)

        if selection:
            array = self.array[selection]
        else:
            array = self.array

        # remove the solvent from the structure if requested
        if del_solvent:
            mask = np.invert(struc.filter_solvent(array))
            if is_stack:
                array = array[:, mask]
            else:
                array = array[mask]

        if del_hydrogen:
            mask = array.element != "H"
            if is_stack:
                array = array[:, mask]
            else:
                array = array[mask]

        obj, frames = _create_object(
            array=array,
            name=name,
            centre=centre,
            style=style,
            collection=collection,
            verbose=verbose,
        )

        if style:
            bl.nodes.create_starting_node_tree(
                object=obj, coll_frames=frames, style=style, color=color
            )

        try:
            obj["entity_ids"] = self.entity_ids
        except AttributeError:
            obj["entity_ids"] = None

        try:
            obj.mn.biological_assemblies = json.dumps(self.assemblies())
        except InvalidFileError:
            obj.mn.biological_assemblies = ""

        if build_assembly and style:
            bl.nodes.assembly_insert(obj)

        # attach the model bpy.Object to the molecule object
        self.object = obj
        # same with the collection of bpy Objects for frames
        self.frames = frames

        return obj

    def assemblies(self, as_array=False):
        """
        Get the biological assemblies of the molecule.

        Parameters
        ----------
        as_array : bool, optional
            Whether to return the assemblies as an array of quaternions. Default is False.

        Returns
        -------
        dict or None
            The biological assemblies of the molecule, as a dictionary of
            transformation matrices, or None if no assemblies are available.
        """
        try:
            assemblies_info = self._assemblies()
        except InvalidFileError:
            return None

        if isinstance(assemblies_info, dict) and as_array:
            return utils.array_quaternions_from_dict(assemblies_info)

        return assemblies_info

    def __repr__(self) -> str:
        """
        Get the string representation of the Molecule object.

        Returns
        -------
        str
            The string representation of the Molecule object.
        """
        return f"<Molecule object: {self.name}>"


def _create_object(
    array,
    name=None,
    centre="",
    style="spherers",
    collection=None,
    world_scale=0.01,
    color_plddt: bool = False,
    verbose=False,
) -> Tuple[bpy.types.Object, bpy.types.Collection]:
    import biotite.structure as struc

    frames = None
    is_stack = isinstance(array, struc.AtomArrayStack)

    try:
        mass = np.array(
            [
                data.elements.get(x, {}).get("standard_mass", 0.0)
                for x in np.char.title(array.element)
            ]
        )
        array.set_annotation("mass", mass)
    except AttributeError as e:
        print(e)

    def centre_array(atom_array, centre):
        if centre == "centroid":
            atom_array.coord -= databpy.centre(atom_array.coord)
        elif centre == "mass":
            atom_array.coord -= databpy.centre(atom_array.coord, weight=atom_array.mass)

    if centre in ["mass", "centroid"]:
        if is_stack:
            for atom_array in array:
                centre_array(atom_array, centre)
        else:
            centre_array(atom_array, centre)

    if is_stack:
        if array.stack_depth() > 1:
            frames = array
        array = array[0]

    if not collection:
        collection = bl.coll.mn()

    bonds_array = []
    bond_idx = []

    if array.bonds:
        bonds_array = array.bonds.as_array()
        bond_idx = bonds_array[:, [0, 1]]
        # the .copy(order = 'C') is to fix a weird ordering issue with the resulting array
        bond_types = bonds_array[:, 2].copy(order="C")

    # creating the blender object and meshes and everything
    bob = databpy.create_bob(
        name=name,
        collection=collection,
        vertices=array.coord * world_scale,
        edges=bond_idx,
    )

    # Add information about the bond types to the model on the edge domain
    # Bond types: 'ANY' = 0, 'SINGLE' = 1, 'DOUBLE' = 2, 'TRIPLE' = 3, 'QUADRUPLE' = 4
    # 'AROMATIC_SINGLE' = 5, 'AROMATIC_DOUBLE' = 6, 'AROMATIC_TRIPLE' = 7
    # https://www.biotite-python.org/apidoc/biotite.structure.BondType.html#biotite.structure.BondType
    if array.bonds:
        bob.store_named_attribute(
            data=bond_types,
            name="bond_type",
            atype=AttributeTypes.INT,
            domain=Domains.EDGE,
        )

    # The attributes for the model are initially defined as single-use functions. This allows
    # for a loop that attempts to add each attibute by calling the function. Only during this
    # loop will the call fail if the attribute isn't accessible, and the warning is reported
    # there rather than setting up a try: except: for each individual attribute which makes
    # some really messy code.

    # I still don't like this as an implementation, and welcome any cleaner approaches that
    # anybody might have.

    def att_atomic_number():
        atomic_number = np.array(
            [
                data.elements.get(x, {"atomic_number": -1}).get("atomic_number")
                for x in np.char.title(array.element)
            ]
        )
        return atomic_number

    def att_atom_id():
        return array.atom_id

    def att_pdb_model_num():
        return array.pdb_model_num

    def att_res_id():
        return array.res_id

    def att_res_name():
        other_res = []
        counter = 0
        id_counter = -1
        res_names = array.res_name
        res_ids = array.res_id
        res_nums = []

        for name in res_names:
            res_num = data.residues.get(name, {"res_name_num": -1}).get("res_name_num")

            if res_num == 9999:
                if (
                    res_names[counter - 1] != name
                    or res_ids[counter] != res_ids[counter - 1]
                ):
                    id_counter += 1

                unique_res_name = str(id_counter + 100) + "_" + str(name)
                other_res.append(unique_res_name)

                num = (
                    np.where(np.isin(np.unique(other_res), unique_res_name))[0][0] + 100
                )
                res_nums.append(num)
            else:
                res_nums.append(res_num)
            counter += 1

        bob.object["ligands"] = np.unique(other_res)
        return np.array(res_nums)

    def att_chain_id():
        if isinstance(array.chain_id[0], int):
            return array.chain_id
        else:
            return np.unique(array.chain_id, return_inverse=True)[1]

    def att_entity_id():
        return array.entity_id

    def att_b_factor():
        return array.b_factor

    def att_occupancy():
        return array.occupancy

    def att_vdw_radii():
        vdw_radii = np.array(
            list(
                map(
                    # divide by 100 to convert from picometres to angstroms which is
                    # what all of coordinates are in
                    lambda x: data.elements.get(x, {}).get("vdw_radii", 100.0) / 100,
                    np.char.title(array.element),
                )
            )
        )
        return vdw_radii * world_scale

    def att_mass():
        return array.mass

    def att_atom_name():
        atom_name = np.array(
            list(map(lambda x: data.atom_names.get(x, -1), array.atom_name))
        )

        return atom_name

    def att_lipophobicity():
        lipo = np.array(
            list(
                map(
                    lambda x, y: data.lipophobicity.get(x, {"0": 0}).get(y, 0),
                    array.res_name,
                    array.atom_name,
                )
            )
        )

        return lipo

    def att_charge():
        charge = np.array(
            list(
                map(
                    lambda x, y: data.atom_charge.get(x, {"0": 0}).get(y, 0),
                    array.res_name,
                    array.atom_name,
                )
            )
        )
        return charge

    def att_color():
        if color_plddt:
            return color.plddt(array.b_factor)
        else:
            return color.color_chains(att_atomic_number(), att_chain_id())

    def att_is_alpha():
        return np.isin(array.atom_name, "CA")

    def att_is_solvent():
        return struc.filter_solvent(array)

    def att_is_backbone():
        """
        Get the atoms that appear in peptide backbone or nucleic acid phosphate backbones.
        Filter differs from the Biotite's `struc.filter_peptide_backbone()` in that this
        includes the peptide backbone oxygen atom, which biotite excludes. Additionally
        this selection also includes all of the atoms from the ribose in nucleic acids,
        and the other phosphate oxygens.
        """
        backbone_atom_names = [
            "N",
            "C",
            "CA",
            "H",  # backbone hydrogen off N
            "HA",  # backbone hydrogen off CA
            "O",  # peptide backbone atoms
            "P",
            "O5'",
            "C5'",
            "C4'",
            "C3'",
            "O3'",  # 'continuous' nucleic backbone atoms
            "O1P",
            "OP1",
            "O2P",
            "OP2",  # alternative names for phosphate O's
            "O4'",
            "C1'",
            "C2'",
            "O2'",  # remaining ribose atoms
        ]

        return np.logical_and(
            np.isin(array.atom_name, backbone_atom_names),
            np.logical_not(struc.filter_solvent(array)),
        )

    def att_is_nucleic() -> npt.NDArray[np.bool_]:
        return struc.filter_nucleotides(array)

    def att_is_peptide() -> npt.NDArray[np.bool_]:
        aa = struc.filter_amino_acids(array)
        con_aa = struc.filter_canonical_amino_acids(array)

        return aa | con_aa

    def att_is_side_chain() -> npt.NDArray[np.bool_]:
        not_backbone = np.logical_not(att_is_backbone())

        return np.logical_and(
            not_backbone, np.logical_or(att_is_nucleic(), att_is_peptide())
        )

    def att_is_hetero():
        return array.hetero

    def att_is_carb() -> npt.NDArray[np.bool_]:
        return struc.filter_carbohydrates(array)

    def att_sec_struct():
        return array.sec_struct

    # these are all of the attributes that will be added to the structure
    # TODO add capcity for selection of particular attributes to include / not include to potentially
    # boost performance, unsure if actually a good idea of not. Need to do some testing.
    attributes = (
        {"name": "res_id", "value": att_res_id, "type": "INT", "domain": "POINT"},
        {"name": "res_name", "value": att_res_name, "type": "INT", "domain": "POINT"},
        {
            "name": "atomic_number",
            "value": att_atomic_number,
            "type": "INT",
            "domain": "POINT",
        },
        {"name": "b_factor", "value": att_b_factor, "type": "FLOAT", "domain": "POINT"},
        {
            "name": "occupancy",
            "value": att_occupancy,
            "type": "FLOAT",
            "domain": "POINT",
        },
        {
            "name": "vdw_radii",
            "value": att_vdw_radii,
            "type": "FLOAT",
            "domain": "POINT",
        },
        {"name": "mass", "value": att_mass, "type": "FLOAT", "domain": "POINT"},
        {"name": "chain_id", "value": att_chain_id, "type": "INT", "domain": "POINT"},
        {
            "name": "pdb_model_num",
            "value": att_pdb_model_num,
            "type": "INT",
            "domain": "POINT",
        },
        {"name": "entity_id", "value": att_entity_id, "type": "INT", "domain": "POINT"},
        {"name": "atom_id", "value": att_atom_id, "type": "INT", "domain": "POINT"},
        {"name": "atom_name", "value": att_atom_name, "type": "INT", "domain": "POINT"},
        {
            "name": "lipophobicity",
            "value": att_lipophobicity,
            "type": "FLOAT",
            "domain": "POINT",
        },
        {"name": "charge", "value": att_charge, "type": "FLOAT", "domain": "POINT"},
        {"name": "Color", "value": att_color, "type": "FLOAT_COLOR", "domain": "POINT"},
        {
            "name": "is_backbone",
            "value": att_is_backbone,
            "type": "BOOLEAN",
            "domain": "POINT",
        },
        {
            "name": "is_side_chain",
            "value": att_is_side_chain,
            "type": "BOOLEAN",
            "domain": "POINT",
        },
        {
            "name": "is_alpha_carbon",
            "value": att_is_alpha,
            "type": "BOOLEAN",
            "domain": "POINT",
        },
        {
            "name": "is_solvent",
            "value": att_is_solvent,
            "type": "BOOLEAN",
            "domain": "POINT",
        },
        {
            "name": "is_nucleic",
            "value": att_is_nucleic,
            "type": "BOOLEAN",
            "domain": "POINT",
        },
        {
            "name": "is_peptide",
            "value": att_is_peptide,
            "type": "BOOLEAN",
            "domain": "POINT",
        },
        {
            "name": "is_hetero",
            "value": att_is_hetero,
            "type": "BOOLEAN",
            "domain": "POINT",
        },
        {"name": "is_carb", "value": att_is_carb, "type": "BOOLEAN", "domain": "POINT"},
        {
            "name": "sec_struct",
            "value": att_sec_struct,
            "type": "INT",
            "domain": "POINT",
        },
    )

    # assign the attributes to the object
    for att in attributes:
        if verbose:
            start = time.process_time()
        try:
            bob.store_named_attribute(
                data=att["value"](),
                name=att["name"],
                atype=att["type"],
                domain=att["domain"],
            )
            if verbose:
                print(f'Added {att["name"]} after {time.process_time() - start} s')
        except Exception as e:
            if verbose:
                print(e)
                warnings.warn(f"Unable to add attribute: {att['name']}")
                print(
                    f'Failed adding {att["name"]} after {time.process_time() - start} s'
                )

    coll_frames = None
    if frames:
        coll_frames = bl.coll.frames(bob.name)
        for i, frame in enumerate(frames):
            frame = databpy.create_object(
                name=bob.name + "_frame_" + str(i),
                collection=coll_frames,
                vertices=frame.coord * world_scale,
            )

    # add custom properties to the actual blender object, such as number of chains, biological assemblies etc
    # currently biological assemblies can be problematic to holding off on doing that
    try:
        bob.object["chain_ids"] = list(np.unique(array.chain_id))
    except AttributeError:
        bob.object["chain_ids"] = None
        warnings.warn("No chain information detected.")

    return bob.object, coll_frames
