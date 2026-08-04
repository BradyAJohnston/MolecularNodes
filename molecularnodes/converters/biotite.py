import MDAnalysis as mda
import numpy as np
from biotite.structure import AtomArray, AtomArrayStack
from MDAnalysis.coordinates.base import SingleFrameReaderBase
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import (
    AtomAttr,
    Atomids,
    Atomnames,
    Bonds,
    ChainIDs,
    Charges,
    Elements,
    ICodes,
    Occupancies,
    Resids,
    Resnames,
    Segids,
    Tempfactors,
)
from MDAnalysis.topology.base import TopologyReaderBase, change_squash


class BiotiteWrapper(object):
    """Biotite Wrapper

    A wrapper for Biotite's AtomArray / AtomArrayStack so that it can be passed
    to MDAnalysis to parse the topology and not consider it as an iterable.

    MDAnalysis supports reading multiple trajectories:
    https://userguide.mdanalysis.org/stable/reading_and_writing.html#reading-multiple-trajectories
    When an iterable is passed as the topology, each of the iterated value is
    parsed for coordinates. Biotite's AtomArray / AtomArrayStack are iterable
    and this leads to a call to parse the individual Atom objects, which is
    incorrect. This wrapper worksaround that issue.

    """

    def __init__(self, structure: AtomArray | AtomArrayStack):
        if not isinstance(structure, (AtomArray, AtomArrayStack)):
            raise ValueError("structure is not an AtomArray or AtomArrayStack")
        self.structure = structure


class BiotiteReader(SingleFrameReaderBase):
    """Biotite Reader

    Read a Biotite AtomArray as a single frame

    """

    format = "BIOTITE"

    units = {"time": None, "length": "Angstrom"}

    @staticmethod
    def _format_hint(thing):
        """Can this reader read *thing*?"""
        return isinstance(thing, BiotiteWrapper)

    def _read_first_frame(self):
        atom_array = self.filename.structure
        if isinstance(atom_array, AtomArrayStack):
            atom_array = atom_array[0]
        self.n_atoms = atom_array.array_length()
        self.ts = ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        ts.positions = atom_array.coord
        ts.frame = 0
        return ts


class BiotiteParser(TopologyReaderBase):
    """Biotite Parser

    Parse a Biotite AtomArray / AtomArrayStack structure to create an MDAnalysis
    Topology. All the mandatory and optional Biotite annotations are converted
    to corresponding MDAnalysis Topology attributes:
    https://biotite-python.org/latest/apidoc/biotite.structure.html

    Only a single (first) AtomArray is used when an AtomArrayStack is passed.

    """

    format = "BIOTITE"

    @staticmethod
    def _format_hint(thing):
        """Can this reader read *thing*?"""
        return isinstance(thing, BiotiteWrapper)

    def parse(self, **kwargs):
        """
        Parse Biotite AtomArray / AtomArrayStack into Topology

        Returns
        -------
        MDAnalysis *Topology* object

        """

        atom_array = self.filename.structure
        if isinstance(atom_array, AtomArrayStack):
            atom_array = atom_array[0]
        n_atoms = atom_array.array_length()

        attrs = []

        # Biotite mandatory annotation categories
        chainids = atom_array.chain_id
        resids = atom_array.res_id
        resnames = atom_array.res_name
        icodes = atom_array.ins_code
        elements = atom_array.element
        # Atom Attr's
        attrs.append(ChainIDs(chainids))
        attrs.append(Atomnames(atom_array.atom_name))
        attrs.append(Elements(elements))
        # Residue Attr's
        residx, (resids, resnames, icodes, chainids) = change_squash(
            (resids, resnames, icodes, chainids),
            (resids, resnames, icodes, chainids),
        )
        n_residues = len(resids)
        attrs.append(Resids(resids))
        attrs.append(Resnames(resnames))
        attrs.append(ICodes(icodes))
        # Segment Attr's
        segidx, (segids,) = change_squash((chainids,), (chainids,))
        n_segments = len(segids)
        attrs.append(Segids(segids))

        # Biotite optional annotation categories
        categories = atom_array.get_annotation_categories()
        for category, Attr in (
            ("atom_id", Atomids),
            ("b_factor", Tempfactors),
            ("occupancy", Occupancies),
            ("charge", Charges),
        ):
            if category not in categories:
                continue
            attrs.append(Attr(atom_array.get_annotation(category)))

        # Carry the connectivity across so that `universe.atoms.bonds` is populated.
        # Biotite stores bonds as (atom_i, atom_j, bond_type); MDAnalysis only wants
        # the index pairs (bond order/type is handled separately on the mesh).
        if atom_array.bonds is not None:
            bond_indices = atom_array.bonds.as_array()[:, :2]
            attrs.append(Bonds([tuple(bond) for bond in bond_indices]))

        return Topology(
            n_atoms,
            n_residues,
            n_segments,
            attrs=attrs,
            atom_resindex=residx,
            residue_segindex=segidx,
        )


class EntityIDs(AtomAttr):
    """Per-atom entity id parsed from the structure file (mmCIF ``label_entity_id``)."""

    attrname = "entity_ids"
    singular = "entity_id"
    dtype = int


class SecStructs(AtomAttr):
    """Per-atom secondary structure code parsed from the structure file."""

    attrname = "sec_structs"
    singular = "sec_struct"
    dtype = int


class Heteros(AtomAttr):
    """Per-atom HETATM flag parsed from the structure file."""

    attrname = "heteros"
    singular = "hetero"
    dtype = bool


class IsCarbs(AtomAttr):
    """Per-atom carbohydrate flag, computed by biotite at parse time.

    Carried across rather than recomputed because it relies on biotite's built-in
    residue knowledge (``filter_carbohydrates``), which has no MDAnalysis equivalent.
    """

    attrname = "is_carbs"
    singular = "is_carb"
    dtype = bool


# Annotations that are read/derived at biotite parse time and cannot be recomputed
# from the topology alone MDAnalysis-side, so they ride across the conversion as custom
# topology attributes. All other (derived) attributes are recomputed by the entity class.
_EXTRA_ANNOTATIONS: dict[str, type[AtomAttr]] = {
    "entity_id": EntityIDs,
    "sec_struct": SecStructs,
    "hetero": Heteros,
    "is_carb": IsCarbs,
}


def universe_from_atoms(
    structure: AtomArray | AtomArrayStack,
) -> mda.Universe:
    """Convert a biotite ``AtomArray``/``AtomArrayStack`` into an MDAnalysis ``Universe``.

    The topology (including bonds) is parsed by :class:`BiotiteParser`. When the input
    is a stack with more than one model, every model is loaded as a trajectory frame so
    multi-model structures (NMR ensembles, multi-model PDB) become multi-frame universes.
    File-parsed annotations that cannot be recomputed from the topology (see
    ``_EXTRA_ANNOTATIONS``) are carried across as custom topology attributes.

    Parameters
    ----------
    structure : AtomArray | AtomArrayStack
        The biotite structure to convert.

    Returns
    -------
    mda.Universe
        A Universe wrapping the structure, with bonds, all coordinate frames, and the
        carried-over file-parsed annotations.
    """
    universe = mda.Universe(BiotiteWrapper(structure))

    # Load every model of a stack as a trajectory frame. `structure.coord` has shape
    # (n_models, n_atoms, 3), exactly what MemoryReader expects.
    if isinstance(structure, AtomArrayStack) and structure.stack_depth() > 1:
        universe.load_new(structure.coord, format=MemoryReader)

    # Carry file-parsed, non-recomputable annotations across as custom attributes.
    reference = structure[0] if isinstance(structure, AtomArrayStack) else structure
    present = reference.get_annotation_categories()
    for name, Attr in _EXTRA_ANNOTATIONS.items():
        if name in present:
            universe.add_TopologyAttr(Attr(reference.get_annotation(name)))

    # Bond order/type is not part of an MDAnalysis topology, and MDAnalysis reorders
    # bonds, so we carry the biotite bond types as an array aligned to the universe's own
    # bond ordering (``atoms.bonds.indices``). The entity stores this on the mesh edge
    # domain when it builds the object. Attached as a plain attribute since it is only
    # needed once, at object-creation time.
    if reference.bonds is not None and hasattr(universe.atoms, "bonds"):
        bond_type_by_pair = {
            frozenset((int(i), int(j))): int(t)
            for i, j, t in reference.bonds.as_array()
        }
        universe._mn_bond_types = np.array(
            [
                bond_type_by_pair.get(frozenset((int(i), int(j))), 0)
                for i, j in universe.atoms.bonds.indices
            ],
            dtype=int,
        )

    return universe
