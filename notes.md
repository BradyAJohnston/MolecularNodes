- pdbx format doesn't have include_bonds as an option (but can be computed afterwards)
- potential for querying the AFDB (maybe leave it to the user to manually download PDB)
- doesn't seem to be able to access assembly transformation matrices
- can't access assemblies from mmtf files

How to get the biological assembly transfrormation matrices (currently via pdbx):

```python
clath = pdbx.PDBxFile.read(rcsb.fetch('1xi4', 'pdbx'))
pdbx.convert._get_transformations(clath.get_category('pdbx_struct_oper_list'))
```

what about doing it with mmtf format??

```python
 clath = mmtf.MMTFFile.read(rcsb.fetch('1xi4', 'mmtf'))

for assembly in clath['bioAssemblyList']:
     for item in assembly:
         if item == 'transformList':
             for matrix in assembly.get(item):
                 print(np.array(matrix.get('matrix')).reshape((4, 4)))

```