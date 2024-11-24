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

        dimensions = np.array(first_line.split(), dtype=int)
        n_atoms = dimensions[0]

        array = np.loadtxt(self.filename, skiprows=1, max_rows=n_atoms, dtype=str)

        attrs = []
        idx = np.arange(1, n_atoms + 1, dtype=int)
        for Attr in (Atomids, Resnums, Resids):
            attrs.append(Attr(idx))
        attrs.append(ChainIDs(array[:, 0].astype(int)))
        attrs.append(Resnames(array[:, 1].astype(str)))

        topo = Topology(n_atoms, n_atoms, 1, attrs=attrs)

        bond_idx = array[:, 2:].astype(int)
        mask = np.logical_and(bond_idx[:, 0] != -1, bond_idx[:, 1] != -1)
        bond_idx = bond_idx[mask, :]

        bonds = Bonds(bond_idx)
        topo.add_TopologyAttr(bonds)

        return topo
