import databpy
import numpy as np
from MDAnalysis import Universe
from MDAnalysis.coordinates.base import ReaderBase
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import (
    Atomids,
    Bonds,
    ChainIDs,
    Resids,
    Resnames,
    # Resnums,
)
from MDAnalysis.lib import util
from MDAnalysis.topology.base import TopologyReaderBase
from ... import color
from ...blender import coll
from ...nodes import nodes
from ..base import EntityType
from .base import Trajectory

DNA_SCALE = 10


class OXDNAParser(TopologyReaderBase):
    def parse(self, **kwargs):
        """Parse topology from an oxDNA topology file.

        Returns
        -------
        top : MDAnalysis.core.topology.Topology
            Topology object
        """
        top = self._parseatoms()

        return top

    @classmethod
    def _is_new_topology(cls, filename) -> bool:
        """Check if topology file is new format.

        Parameters
        ----------
        filename : str
            Path to topology file

        Returns
        -------
        bool
            True if topology file is new format
        """
        with open(filename) as f:
            return "5->3" in f.readline()

    @classmethod
    def _read_topo_new(cls, filename) -> Topology:
        """Read topology from new format oxDNA topology file.

        Parameters
        ----------
        filename : str
            Path to topology file

        Returns
        -------
        topo : MDAnalysis.core.topology.Topology
            Topology object
        """
        with open(filename) as f:
            lines = f.readlines()
        n_atoms, n_chains, direction = np.array(lines[0].split())
        n_atoms = int(n_atoms)
        atom_idx = np.arange(n_atoms)

        res_name_list = []
        chain_id_list = []
        for i, line in enumerate(lines[1:]):
            is_rna = "type=RNA" in line
            _is_dna = not is_rna

            line_split = line.split()
            bases = line_split[0]

            base_list = []
            start = 0
            end = 0
            in_custom_base = False
            for i, letter in enumerate(bases):
                if letter == "(":
                    start = i + 1
                    in_custom_base = True
                    continue
                if letter == ")":
                    end = i
                    base_list.append(line[start:end])
                    in_custom_base = False
                    continue

                if in_custom_base:
                    continue

                base_list.append(letter)
            chain_id_list.append(np.repeat(i, len(base_list)))

            res_name_list.append(np.array(base_list))

        res_names = np.hstack(res_name_list)
        chain_ids = np.hstack(chain_id_list)

        bond_idx = np.zeros((n_atoms, 2), dtype=int)
        bond_idx[:, :] = -1

        for i in atom_idx:
            if i == 1:
                continue
            if chain_ids[i] == chain_ids[-1]:
                bond_idx[i, :] = np.array((i, i - 1), dtype=int)

        mask = np.logical_and(bond_idx[:, 0] != -1, bond_idx[:, 1] != -1)
        bond_idx = bond_idx[mask, :]

        attrs = [
            Atomids(atom_idx),
            Resids(atom_idx + 1),
            ChainIDs(chain_ids),
            Resnames(res_names),
        ]

        topo = Topology(
            n_atoms=n_atoms, n_res=n_atoms, attrs=attrs, atom_resindex=atom_idx
        )
        topo.add_TopologyAttr(Bonds(bond_idx))

        return topo

    @classmethod
    def _read_topo_old(cls, filename) -> Topology:
        """Read topology from old format oxDNA topology file.

        Parameters
        ----------
        filename : str
            Path to topology file

        Returns
        -------
        topo : MDAnalysis.core.topology.Topology
            Topology object
        """
        with open(filename) as f:
            first_line = f.readline()

        dimensions = np.array(first_line.split())
        n_atoms = int(dimensions[0])
        array = np.loadtxt(filename, skiprows=1, max_rows=n_atoms, dtype=str)

        # each topology item has two bond columns, which say what the base is bonded
        # from and what it is bonded to. -1 means it is not bonded
        # we need to turn that into an array of bond pairs
        bond_idx = np.zeros((n_atoms * 2, 2), int)
        row_numbers = np.arange(n_atoms, dtype=int)
        for i in range(2):
            rows = row_numbers + (i * n_atoms)
            bond_idx[rows, 0] = array[:, (i + 2)]
            bond_idx[rows, 1] = row_numbers

        # drop any that have -1 as they aren't bonded to anything
        mask = np.logical_and(bond_idx[:, 0] != -1, bond_idx[:, 1] != -1)
        bond_idx = bond_idx[mask, :]

        atom_idx = np.arange(n_atoms, dtype=int)
        attrs = [
            Atomids(atom_idx),
            Resids(atom_idx + 1),
            ChainIDs(array[:, 0].astype(int)),
            Resnames(array[:, 1].astype(str)),
        ]

        topo = Topology(
            n_atoms=n_atoms, n_res=n_atoms, attrs=attrs, atom_resindex=atom_idx
        )
        topo.add_TopologyAttr(Bonds(bond_idx))

        return topo

    def _parseatoms(self):
        """Parse atoms from topology file.

        Returns
        -------
        topo : MDAnalysis.core.topology.Topology
            Topology object
        """
        if self._is_new_topology(self.filename):
            return self._read_topo_new(self.filename)
        else:
            return self._read_topo_old(self.filename)


def _is_info_line(line: str):
    """
    Check if a line contains simulation information.

    Parameters
    ----------
    line : str
        Line from the OXDNA trajectory file.

    Returns
    -------
    bool
        True if line contains simulation information, False otherwise.
    """
    return line.startswith("t = ") or line.startswith("b = ") or line.startswith("E = ")


class OXDNAReader(ReaderBase):
    """
    Reader for OXDNA trajectory files.

    Parameters
    ----------
    filename : str
        Path to the OXDNA trajectory file.
    **kwargs : dict
        Additional arguments to pass to the reader, must include 'n_atoms'.

    Attributes
    ----------
    n_atoms : int
        Number of atoms in the system.
    ts : Timestep
        Current timestep object.
    n_frames : int
        Total number of frames in trajectory.
    """

    def __init__(self, filename, **kwargs):
        super(OXDNAReader, self).__init__(filename, **kwargs)

        self.n_atoms = kwargs["n_atoms"]
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)

        oxdnafile = self._oxdnafile = util.anyopen(filename, "rb")

        starts = []
        stops = []

        line = "a new line"
        previous_is_info = False
        current_is_info = False
        while line:
            line = oxdnafile.readline().decode()
            current_is_info = _is_info_line(line)

            # get the lines where the data values start and stop, which
            # are separated by the several lines of overall information
            if not current_is_info and previous_is_info:
                starts.append(oxdnafile.tell() - len(line))
            if current_is_info and not previous_is_info:
                stops.append(oxdnafile.tell() - len(line))

            previous_is_info = current_is_info

        # add the end of the file
        stops += [oxdnafile.tell()]

        self._start_offsets = starts
        self._stop_offsets = stops[1:]  # drop the first
        self.n_frames = len(self._start_offsets)

        self._read_frame(0)

    def _reopen(self):
        """
        Reopen the trajectory file and reset the timestep frame.
        """
        self.close()
        self._oxdnafile = util.anyopen(self.filename, "rb")
        self.ts.frame = -1

    def _read_next_timestep(self, ts=None):
        """
        Read the next timestep from the trajectory.

        Parameters
        ----------
        ts : Timestep, optional
            If provided, update this timestep object instead of creating a new one.

        Returns
        -------
        Timestep
            The timestep object containing the frame data.
        """
        frame = self.frame + 1
        return self._read_frame(frame)

    def _read_frame(self, frame):
        """
        Read a specific frame from the trajectory.

        Parameters
        ----------
        frame : int
            Frame number to read.

        Returns
        -------
        Timestep
            The timestep object containing the frame data.

        Raises
        ------
        OSError
            If the frame number is out of range.
        """
        try:
            start = self._start_offsets[frame]
            stop = self._stop_offsets[frame]
        except IndexError:
            raise OSError from None

        self._oxdnafile.seek(start)
        chunk = self._oxdnafile.read(stop - start)

        array = np.array(
            [
                np.array(line.split(), dtype=float)
                for line in chunk.decode().splitlines()
            ]
        )

        # TODO: also access and update the other values
        self.ts.positions = array[:, :3]

        for i, name in enumerate(
            ("base_vector", "base_normal", "velocity", "angular_velocity")
        ):
            starting_column = 3 * (i + 1)
            if starting_column >= array.shape[1]:
                continue
            cols = np.arange(3) + starting_column
            self.ts.data[name] = array[:, cols]
        self.ts.frame = frame

        return self.ts

    def close(self):
        """
        Close the trajectory file.
        """
        self._oxdnafile.close()


class OXDNA(Trajectory):
    """
    A class to handle oxDNA trajectory data.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        The MDAnalysis Universe object containing the trajectory data
    world_scale : float, optional
        Scaling factor for the world coordinates, by default 0.01

    Attributes
    ----------
    _entity_type : EntityType
        Type of the molecular entity
    _att_names : tuple
        Names of the attributes to track
    """

    def __init__(self, universe: Universe, world_scale: float = 0.01):
        super().__init__(universe=universe, world_scale=world_scale * DNA_SCALE)
        self._entity_type = EntityType.MD_OXDNA
        self._att_names = (
            "base_vector",
            "base_normal",
            "velocity",
            "angular_velocity",
        )

    def _create_object(
        self, style: str | None = "oxdna", name: str = "NewUniverseObject"
    ):
        """
        Create a new object with the trajectory data.

        Parameters
        ----------
        style : str, optional
            Style of the object representation, by default "oxdna"
        name : str, optional
            Name of the new object, by default "NewUniverseObject"

        Returns
        -------
        bpy.types.Object
            The created Blender object
        """
        self.object = databpy.create_object(
            name=name,
            collection=coll.mn(),
            vertices=self.univ_positions,
            edges=self.bonds,
        )
        self._update_timestep_values()

        for name in ("chain_id", "res_id", "res_name"):
            if name == "res_name":
                att_name = "res_num"
            else:
                att_name = name
            self.store_named_attribute(getattr(self, att_name), name)

        self.store_named_attribute(
            data=color.color_chains_equidistant(self.chain_id),
            name="Color",
            atype=databpy.AttributeTypes.FLOAT_COLOR,
        )

        if style:
            nodes.create_starting_node_tree(self.object, style="oxdna", color=None)

        return self.object

    def set_frame(self, frame: int) -> None:
        super()._update_positions(frame)
        self._update_timestep_values()

    def _update_timestep_values(self):
        """
        Update the timestep values for all tracked attributes.
        """
        for name in self._att_names:
            try:
                self.store_named_attribute(
                    self.universe.trajectory.ts.data[name] * self.world_scale, name=name
                )
            except KeyError as e:
                print(e)
