import itertools
import re
import bpy
import databpy as db
import numpy as np
import pytest
import molecularnodes as mn
from molecularnodes.entities.base import EntityType
from molecularnodes.entities.molecule.pdbx import MultipleStructuresError
from molecularnodes.entities.reload import reload_entity
from molecularnodes.nodes.geometry import (
    StyleBallAndStick,
    StyleCartoon,
    StyleRibbon,
    StyleSpheres,
    StyleSurface,
)
from molecularnodes.session import get_session
from .constants import codes, data_dir
from .utils import GeometrySet, NumpySnapshotExtension, mn_ops

STYLES_TO_TEST = [
    "cartoon",
    "ribbon",
    "spheres",
    "surface",
    "ball_and_stick",
]

# `surface` builds its mesh through a volume, and the marching cubes step lands on a
# slightly different vertex count on each platform's Blender build. It is also by far
# the slowest style, so it is only exercised on the smallest structure.
SURFACE_CODE = "1BNA"
STYLE_PARAMS = [
    (code, assembly, style)
    for code, assembly, style in itertools.product(codes, [True, False], STYLES_TO_TEST)
    if style != "surface" or code == SURFACE_CODE
]


@pytest.mark.parametrize("code, assembly, style", STYLE_PARAMS)
def test_style_1(fetch, snapshot, code, assembly, style):
    mol = fetch(code)
    with mol.tree.reset() as (atoms, join):
        match style:
            case "ball_and_stick":
                style_node = StyleBallAndStick(sphere="Mesh")
            case "spheres":
                style_node = StyleSpheres(sphere="Mesh")
            case "cartoon":
                style_node = StyleCartoon()
            case "ribbon":
                style_node = StyleRibbon()
            case "surface":
                style_node = StyleSurface()

        assembly = (
            mn.nodes.geometry.AssemblyInstance(data_object=mol.create_data_object())
            if assembly
            else None
        )
        (atoms >> style_node >> assembly >> join)

    assert snapshot == GeometrySet(mol.object).summary()


@pytest.mark.parametrize(
    "code, format", list(itertools.product(codes, ["bcif", "cif", "pdb"]))
)
def test_download_format(fetch, code, format):
    mol = fetch(code, format=format)
    assert mol.props.entity_type == EntityType.MOLECULE.value
    with db.ObjectTracker() as o:
        mn_ops.import_molecule(code=code, file_format=format, cache_dir=str(data_dir))
        mol2 = get_session().match(o.latest())
    assert isinstance(mol2, mn.Molecule)

    assert np.allclose(mol.position, mol2.position)


@pytest.mark.parametrize("code", codes)
def test_style_positions(fetch, snapshot_custom: NumpySnapshotExtension, code):
    mol = fetch(code)
    assert snapshot_custom == mol.position


def test_local_pdb(fetch, snapshot_custom):
    molecules = [mn.Molecule.load(data_dir / f"1l58.{ext}") for ext in ("cif", "pdb")]
    molecules.append(fetch("1l58", format="bcif"))
    for mol in molecules:
        assert snapshot_custom == mol.named_attribute("position")


def test_pdb_blank_res_id(snapshot):
    # VESTA writes PDB files with blank res_id columns, which biotite refuses
    # to parse without help (#1091)
    mol = mn.Molecule.load(data_dir / "vesta_blank_res_id.pdb")
    assert mol.atoms.n_atoms == 20
    assert (mol.atoms.resids == 1).all()
    assert snapshot == mol.position


def test_cif_no_group_pdb(snapshot):
    # gemmi omits the optional group_PDB column from atom_site, which biotite
    # requires for the hetero annotation (#19)
    mol = mn.Molecule.load(data_dir / "1BNA_gemmi.cif")
    assert mol.atoms.n_atoms == 566
    assert snapshot == mol.position


@pytest.mark.filterwarnings("ignore:.*elements were guessed.*:UserWarning")
def test_pdb_no_bonds(snapshot):
    mol = mn.Molecule.load(data_dir / "no_bonds.pdb")
    assert isinstance(mol.object.data, bpy.types.Mesh)
    assert len(mol.object.data.edges) == 0
    assert snapshot == mol.position


def test_rcsb_nmr(fetch, snapshot_custom):
    mol = fetch("2M6Q")
    # multi-model (NMR) structures now load as trajectory frames of the Universe
    assert mol.universe.trajectory.n_frames > 1

    mol.add_style("cartoon")
    assert snapshot_custom == mol.named_attribute("position")

    # advancing the frame pushes the corresponding model's positions onto the mesh
    mol.set_frame(0)
    pos_0 = mol.named_attribute("position")
    mol.set_frame(5)
    pos_5 = mol.named_attribute("position")
    assert not np.allclose(pos_0, pos_5)


def test_load_small_mol(snapshot):
    mol = mn.Molecule.load(data_dir / "ASN.cif")
    assert mol.props.entity_type == EntityType.MOLECULE.value
    assert snapshot == GeometrySet(mol.object, strict=True).summary()


def _terminal_backbone_indices(mol) -> tuple[int, int]:
    """Indices of the first protein residue's backbone N and the last one's C."""
    protein = mol.universe.select_atoms("protein")
    first_n = protein[(protein.resids == protein.resids.min()) & (protein.names == "N")]
    last_c = protein[(protein.resids == protein.resids.max()) & (protein.names == "C")]
    return int(first_n[0].ix), int(last_c[0].ix)


def test_cyclic_peptide_closing_bond(fetch):
    # 1KAL is a head-to-tail cyclic peptide. The closing C -> N bond is absent
    # from the file's struct_conn (which lists only the disulfides), so it is
    # detected geometrically on import and must appear in the bonds and the
    # mesh edges (#858)
    mol = fetch("1KAL")
    first_n, last_c = _terminal_backbone_indices(mol)

    bonds = np.sort(mol.atoms.bonds.indices, axis=1)
    assert np.any(np.all(bonds == sorted((first_n, last_c)), axis=1))

    edges = np.sort(np.array([edge.vertices for edge in mol.object.data.edges]), axis=1)
    assert np.any(np.all(edges == sorted((first_n, last_c)), axis=1))

    # a linear peptide must not gain a closing bond
    mol_linear = fetch("4ozs")
    first_n, last_c = _terminal_backbone_indices(mol_linear)
    bonds = np.sort(mol_linear.atoms.bonds.indices, axis=1)
    assert not np.any(np.all(bonds == sorted((first_n, last_c)), axis=1))


RESTRAINT_BLOCK = """data_link_SS
loop_
  _chem_link_bond.link_id
  _chem_link_bond.atom_1_comp_id
  _chem_link_bond.atom_id_1
  _chem_link_bond.atom_2_comp_id
  _chem_link_bond.atom_id_2
  _chem_link_bond.type
  _chem_link_bond.value_dist
  _chem_link_bond.value_dist_esd
   SS 1 SG 2 SG single 2.031 0.020
"""


def test_cif_phenix_restraint_block(tmp_path):
    # phenix.real_space_refine indents the tags under `loop_` and appends
    # restraint definition blocks with no atoms after the model
    model = re.sub(
        r"^(_|loop_)",
        r"  \1",
        (data_dir / "1BNA.cif").read_text(),
        flags=re.MULTILINE,
    )
    path = tmp_path / "refined.cif"
    path.write_text(model.replace("  loop_", "loop_") + "\n" + RESTRAINT_BLOCK)
    mol = mn.Molecule.load(path)
    single = mn.Molecule.load(data_dir / "1BNA.cif")
    assert mol.atoms.n_atoms == single.atoms.n_atoms
    assert np.allclose(mol.position, single.position)


def test_cif_multiple_structure_blocks(tmp_path):
    path = tmp_path / "two.cif"
    path.write_text(
        (data_dir / "1BNA.cif").read_text() + "\n" + (data_dir / "1l58.cif").read_text()
    )
    with pytest.raises(MultipleStructuresError) as e:
        mn.Molecule.load(path)
    assert e.value.blocks == ["1BNA", "1L58"]

    for block, code in zip(e.value.blocks, ["1BNA.cif", "1l58.cif"]):
        mol = mn.Molecule.load(path, data_block=block)
        single = mn.Molecule.load(data_dir / code)
        assert mol.atoms.n_atoms == single.atoms.n_atoms
        assert mol.props.data_block == block
        # reloading reads the same block back
        reloaded = reload_entity(mol.object)
        assert isinstance(reloaded, mn.Molecule)
        assert reloaded.atoms.n_atoms == single.atoms.n_atoms
