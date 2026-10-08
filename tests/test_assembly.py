import itertools
from os.path import dirname, join, realpath
import biotite.structure.io.pdb as biotite_pdb
import biotite.structure.io.pdbx as biotite_cif
import numpy as np
import pytest
import molecularnodes as mn
import molecularnodes.entities.molecule.pdb as pdb
import molecularnodes.entities.molecule.pdbx as pdbx
from molecularnodes.download import StructureDownloader

DATA_DIR = join(dirname(realpath(__file__)), "data")


@pytest.fixture(scope="module")
def path_4ins():
    return StructureDownloader().download("4INS", format="pdb")


@pytest.mark.parametrize(
    "pdb_id, format", list(itertools.product(["1f2n", "5zng"], ["pdb", "cif"]))
)
def test_get_transformations(pdb_id, format):
    """
    Compare an assembly built from transformation information in
    MolecularNodes with assemblies built in Biotite.
    """
    path = join(DATA_DIR, f"{pdb_id}.{format}")
    if format == "pdb":
        pdb_file = biotite_pdb.PDBFile.read(path)
        atoms = biotite_pdb.get_structure(pdb_file, model=1)
        ref_assembly = biotite_pdb.get_assembly(pdb_file, model=1)
        test_parser = pdb.PDBAssemblyParser(pdb_file)
    elif format == "cif":
        cif_file = biotite_cif.CIFFile().read(path)
        atoms = biotite_cif.get_structure(
            # Make sure `label_asym_id` is used instead of `auth_asym_id`
            cif_file,
            model=1,
            use_author_fields=False,
        )
        ref_assembly = biotite_cif.get_assembly(cif_file, model=1)
        test_parser = pdbx.CIFAssemblyParser(cif_file)
    else:
        raise ValueError(f"Format '{format}' does not exist")

    assembly_id = test_parser.list_assemblies()[0]
    test_transformations = test_parser.get_transformations(assembly_id)

    check_transformations(test_transformations, atoms, ref_assembly)


@pytest.mark.parametrize("assembly_id", [str(i + 1) for i in range(5)])
def test_get_transformations_cif(assembly_id):
    """
    Compare an assembly built from transformation information in
    MolecularNodes with assemblies built in Biotite.

    In this case all assemblies from a structure with more complex
    operation expressions are tested
    """
    cif_file = biotite_cif.CIFFile().read(join(DATA_DIR, "1f2n.cif"))
    atoms = biotite_cif.get_structure(
        # Make sure `label_asym_id` is used instead of `auth_asym_id`
        cif_file,
        model=1,
        use_author_fields=False,
    )
    ref_assembly = biotite_cif.get_assembly(cif_file, model=1, assembly_id=assembly_id)

    test_parser = pdbx.CIFAssemblyParser(cif_file)

    test_transformations = test_parser.get_transformations(assembly_id)

    check_transformations(test_transformations, atoms, ref_assembly)


@pytest.mark.parametrize("assembly_id", [str(i + 1) for i in range(7)])
def test_get_transformations_pdb(assembly_id, path_4ins):
    """
    Compare each assembly of a PDB file containing multiple BIOMOLECULE
    blocks with non-identity rotations against the assemblies built in Biotite.
    """
    pdb_file = biotite_pdb.PDBFile.read(str(path_4ins))
    atoms = biotite_pdb.get_structure(pdb_file, model=1)
    ref_assembly = biotite_pdb.get_assembly(pdb_file, model=1, assembly_id=assembly_id)

    test_parser = pdb.PDBAssemblyParser(pdb_file)
    assert test_parser.list_assemblies() == [str(i + 1) for i in range(7)]

    test_transformations = test_parser.get_transformations(assembly_id)

    check_transformations(test_transformations, atoms, ref_assembly)


@pytest.mark.parametrize("format", ["pdb", "cif"])
def test_assemblies_as_array(format):
    """Assemblies from both formats convert to the per-chain transform array."""
    mol = mn.Molecule.load(join(DATA_DIR, f"1cd3.{format}"))
    transforms = mol.assemblies(as_array=True)
    assert transforms is not None
    assert set(np.unique(transforms["assembly_id"])) == set(
        range(1, len(mol.assemblies()) + 1)
    )


def test_assemblies_parse_error_warns(caplog, tmp_path):
    """A malformed REMARK 350 is reported, rather than looking like no assemblies."""
    with open(join(DATA_DIR, "1f2n.pdb")) as f:
        lines = f.read().splitlines()
    # drop a single BIOMT line so the transformation vectors no longer come in threes
    biomt = [
        i for i, line in enumerate(lines) if line[11:].lstrip().startswith("BIOMT")
    ]
    del lines[biomt[0]]
    path = tmp_path / "1f2n.pdb"
    path.write_text("\n".join(lines))
    reader = pdb.PDBReader(path)

    with caplog.at_level("WARNING"):
        assert reader.assemblies() == ""
    assert "Failed to parse biological assemblies" in caplog.text


def test_assemblies_none_is_silent(caplog, tmp_path):
    """A PDB file without assembly records returns no assemblies and no warning."""
    with open(join(DATA_DIR, "1f2n.pdb")) as f:
        lines = [
            line
            for line in f.read().splitlines()
            if not line.startswith(("REMARK 300", "REMARK 350"))
        ]
    path = tmp_path / "1f2n.pdb"
    path.write_text("\n".join(lines))
    reader = pdb.PDBReader(path)

    with caplog.at_level("WARNING"):
        assert reader.assemblies() == {}
    assert "biological assemblies" not in caplog.text


def check_transformations(transformations, atoms, ref_assembly):
    """
    Check if the given transformations applied on the given atoms
    results in the given reference assembly.
    """
    test_assembly = None
    for transformation in transformations:
        chain_ids = transformation["chain_ids"]
        matrix = np.array(transformation["matrix"])
        translation = matrix[:3, 3]
        rotation = matrix[:3, :3]
        sub_assembly = atoms[np.isin(atoms.chain_id, chain_ids)].copy()
        sub_assembly.coord = np.dot(rotation, sub_assembly.coord.T).T
        sub_assembly.coord += translation
        if test_assembly is None:
            test_assembly = sub_assembly
        else:
            test_assembly += sub_assembly

    assert test_assembly.array_length() == ref_assembly.array_length()
    # The atom name is used as indicator of correct atom ordering here
    assert np.all(test_assembly.atom_name == ref_assembly.atom_name)
    assert np.allclose(test_assembly.coord, ref_assembly.coord, atol=1e-4)
