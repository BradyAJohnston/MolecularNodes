from biotite.structure import AtomArray, AtomArrayStack
from MDAnalysis.coordinates.base import SingleFrameReaderBase
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import (
    Atomids,
    Atomnames,
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

        return Topology(
            n_atoms,
            n_residues,
            n_segments,
            attrs=attrs,
            atom_resindex=residx,
            residue_segindex=segidx,
        )
