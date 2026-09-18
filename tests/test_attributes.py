import itertools
import numpy as np
import pytest
from MDAnalysis.tests.datafiles import DCD, GRO, PSF, XTC
import molecularnodes as mn
from molecularnodes.nodes.geometry import StyleRibbon, StyleSpheres
from .constants import codes, data_dir
from .utils import GeometrySet

formats = ["pdb", "cif", "bcif"]


@pytest.mark.parametrize("code, format", list(itertools.product(codes, formats)))
def test_attribute(snapshot, code, format):
    mol = mn.Molecule.fetch(code, cache=data_dir, format=format)
    # no style applied, so the geometry is the parsed file and is exact on all platforms
    assert snapshot == GeometrySet(mol.object, strict=True).summary()


def test_store_named_attribute(snapshot_custom):
    mol = mn.Molecule.fetch("8H1B", cache=data_dir, format="bcif")
    before = mol.named_attribute("position")
    mol.store_named_attribute(mol.named_attribute("position") + 10, "position")
    after = mol.named_attribute("position")

    assert not np.allclose(before, after)


def test_uv_map(snapshot_custom):
    mol = mn.Molecule.fetch("1cd3", cache=data_dir, format="bcif")
    with mol.tree.reset() as (atoms, join):
        atoms >> StyleRibbon(uv_map=True, quality=1) >> join
    assert snapshot_custom == mol.named_attribute("uv_map", evaluate=True)[:1000]
    assert snapshot_custom == mol.named_attribute("uv_map", evaluate=True)[-1000:]


def test_bond_attributes(snapshot):
    mol = mn.Molecule.fetch("1BNA", cache=data_dir, format="bcif")
    with mol.tree.reset() as (atoms, join):
        atoms >> StyleSpheres(sphere="Mesh") >> join

    assert snapshot == GeometrySet(mol.object).summary()


@pytest.mark.parametrize("format", formats)
def test_charge_from_file(format):
    # 8U8W carries formal charges on its ions (NA 1+, 2x CL 1-, 2x IOD 1-) in the PDB
    # charge column and in `pdbx_formal_charge`; they are stored verbatim, with 0 on
    # every other atom
    mol = mn.Molecule.fetch("8U8W", cache=data_dir, format=format)
    charge = mol.named_attribute("charge")
    is_ion = np.isin(mol.universe.atoms.resnames, ["NA", "CL", "IOD"])
    assert is_ion.sum() == 5
    assert np.array_equal(np.sort(charge[is_ion]), [-1, -1, -1, -1, 1])
    assert np.all(charge[~is_ion] == 0)


@pytest.mark.parametrize("format", formats)
def test_charge_absent_without_file_charges(format):
    # 4ozs has no charges in the file (blank column / all `?`), so no `charge`
    # attribute is stored rather than a looked-up or zero-filled one
    mol = mn.Molecule.fetch("4ozs", cache=data_dir, format=format)
    assert "charge" not in mol.list_attributes()


def test_charge_pdb_column(tmp_path):
    # the PDB charge column (cols 79-80) is stored verbatim, blank entries as 0
    atoms = [("N", "N", "1+"), ("CA", "C", "  "), ("C", "C", "  "), ("O", "O", "1-")]
    lines = [
        f"ATOM  {i + 1:5d}  {name:<3s} ALA A   1    {i:8.3f}{0:8.3f}{0:8.3f}  1.00 10.00          {element:>2s}{charge}"
        for i, (name, element, charge) in enumerate(atoms)
    ]
    path = tmp_path / "charged.pdb"
    path.write_text("\n".join(lines) + "\nEND\n")
    mol = mn.Molecule.load(path)
    assert np.array_equal(mol.named_attribute("charge"), [1.0, 0.0, 0.0, -1.0])


def test_charge_topology():
    # a simulation topology with partial charges keeps them
    mol = mn.Molecule.load(PSF, DCD)
    charge = mol.named_attribute("charge")
    assert np.any(charge != 0)
    assert np.allclose(charge, mol.universe.atoms.charges)


def test_charge_absent_without_topology_charges():
    # a topology without charges stores no `charge` attribute rather than zeros
    mol = mn.Molecule.load(GRO, XTC)
    assert not hasattr(mol.universe.atoms, "charges")
    assert "charge" not in mol.list_attributes()
