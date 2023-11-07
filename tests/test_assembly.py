from os.path import join, dirname, realpath
import itertools
import pytest
import numpy as np
import biotite.structure.io.pdb as biotite_pdb
import biotite.structure.io.pdbx as biotite_cif
import biotite.structure.io.mmtf as biotite_mmtf
import molecularnodes.assembly.pdb as pdb
import molecularnodes.assembly.cif as cif
import molecularnodes.assembly.mmtf as mmtf


DATA_DIR = join(dirname(realpath(__file__)), "data")


@pytest.mark.parametrize("pdb_id, file_format", itertools.product(
    ["1f2n", "5zng"],
    ["pdb", "cif", "mmtf"]
))
def test_get_transformations(pdb_id, file_format):
    """
    Compare an assembly built from transformation information in
    MolecularNodes with assemblies built in Biotite.
    """
    path = join(DATA_DIR, f"{pdb_id}.{file_format}")
    if file_format == "pdb":
        pdb_file = biotite_pdb.PDBFile.read(path)
        atoms = biotite_pdb.get_structure(pdb_file, model=1)
        ref_assembly = biotite_pdb.get_assembly(pdb_file, model=1)
        test_parser = pdb.PDBAssemblyParser(pdb_file)
    elif file_format == "cif":
        cif_file = biotite_cif.PDBxFile.read(path)
        atoms = biotite_cif.get_structure(
            # Make sure `label_asym_id` is used instead of `auth_asym_id`
            cif_file, model=1, use_author_fields=False
        )
        ref_assembly = biotite_cif.get_assembly(cif_file, model=1)
        test_parser = cif.CIFAssemblyParser(cif_file)
    elif file_format == "mmtf":
        mmtf_file = biotite_mmtf.MMTFFile.read(path)
        atoms = biotite_mmtf.get_structure(mmtf_file, model=1)
        try:
            ref_assembly = biotite_mmtf.get_assembly(mmtf_file, model=1)
        except NotImplementedError:
            pytest.skip(
                "The limitation of the function does not support this "
                "structure"
            )
        test_parser = mmtf.MMTFAssemblyParser(mmtf_file) 
    else:
        raise ValueError(f"Format '{file_format}' does not exist")
    
    assembly_id = test_parser.list_assemblies()[0]
    test_transformations = test_parser.get_transformations(assembly_id)
    
    check_transformations(test_transformations, atoms, ref_assembly)


@pytest.mark.parametrize("assembly_id", [str(i+1) for i in range(6)])
def test_get_transformations_cif(assembly_id):
    """
    Compare an assembly built from transformation information in
    MolecularNodes with assemblies built in Biotite.

    In this case all assemblies from a structure with more complex
    operation expressions are tested
    """
    cif_file = biotite_cif.PDBxFile.read(join(DATA_DIR, "1f2n.cif"))
    atoms = biotite_cif.get_structure(
        # Make sure `label_asym_id` is used instead of `auth_asym_id`
        cif_file, model=1, use_author_fields=False
    )
    ref_assembly = biotite_cif.get_assembly(
        cif_file, model=1, assembly_id=assembly_id
    )
    
    test_parser = cif.CIFAssemblyParser(cif_file)
    test_transformations = test_parser.get_transformations(assembly_id)
    
    check_transformations(test_transformations, atoms, ref_assembly)


def check_transformations(transformations, atoms, ref_assembly):
    """
    Check if the given transformations applied on the given atoms
    results in the given reference assembly.
    """
    test_assembly = None
    for chain_ids, rotation, translation in transformations:
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