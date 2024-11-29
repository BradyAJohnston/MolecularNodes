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

    @classmethod
    def _is_new_topology(self, filename) -> bool:
        with open(filename) as f:
            return "5->3" in f.readline()

    @classmethod
    def _read_topo_new(cls, filename) -> Topology:
        with open(filename) as f:
            lines = f.readlines()
        n_atoms, n_chains, direction = np.array(lines[0].split())
        n_atoms = int(n_atoms)
        atom_idx = np.arange(n_atoms)

        res_name_list = []
        chain_id_list = []
        for i, line in enumerate(lines[1:]):
            is_rna = "type=RNA" in line
            is_dna = not is_rna

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
        if self._is_new_topology(self.filename):
            return self._read_topo_new(self.filename)
        else:
            return self._read_topo_old(self.filename)
