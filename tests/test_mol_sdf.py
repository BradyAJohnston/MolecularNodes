def open_structure_local_mol(file_path):
    from biotite.structure.io.mol import MOLFile

    file = MOLFile.read(file_path)

    mol = file.get_structure()

    return mol,file

file_path = 'data/caffeine-3D-structure-CT1001987571.sdf'

mol,file = open_structure_local_mol(file_path)

file_path = 'data/caffeine.mol'

mol,file = open_structure_local_mol(file_path)

print('success')
