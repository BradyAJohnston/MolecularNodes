import numpy as np
from MDAnalysis.topology.base import TopologyReaderBase
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import (
    Atomids,
    Bonds,
    ChainIDs,
    Resids,
    Resnames,
    Resnums,
)


class OXDNAParser(TopologyReaderBase):
    def parse(self, **kwargs):
        top = self._parseatoms()

        return top

    def _parseatoms(self):
        with open(self.filename) as f:
            first_line = f.readline()

        is_new_topology = "5->3" in first_line

        if not is_new_topology:
            dimensions = np.array(first_line.split(), dtype=int)
            n_atoms = dimensions[0]
            array = np.loadtxt(self.filename, skiprows=1, max_rows=n_atoms, dtype=str)

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

        attrs = []
        idx = np.arange(1, n_atoms + 1, dtype=int)
        for Attr in (Atomids, Resnums, Resids):
            attrs.append(Attr(idx))
        attrs.append(ChainIDs(array[:, 0].astype(int)))
        attrs.append(Resnames(array[:, 1].astype(str)))

        topo = Topology(n_atoms, n_atoms, 1, attrs=attrs)

        bonds = Bonds(bond_idx)
        topo.add_TopologyAttr(bonds)

        return topo
