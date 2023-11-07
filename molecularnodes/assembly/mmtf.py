import numpy as np
from . import AssemblyParser

class MMTFAssemblyParser(AssemblyParser):
    ### Implementation adapted from ``biotite.structure.io.mmtf.assembly``

    def __init__(self, mmtf_file):
        self._file = mmtf_file
    

    def list_assemblies(self):
        import biotite.structure.io.mmtf as mmtf
        return mmtf.list_assemblies(self._file)
    

    def get_transformations(self, assembly_id):
        import biotite
        # Find desired assembly
        selected_assembly = None
        if not "bioAssemblyList" in self._file:
            raise biotite.InvalidFileError(
                "File does not contain assembly information "
                "(missing 'bioAssemblyList')"
            )
        for assembly in self._file["bioAssemblyList"]:
            current_assembly_id = assembly["name"]
            transform_list = assembly["transformList"]
            if current_assembly_id == assembly_id:
                selected_assembly = transform_list
                break
        if selected_assembly is None:
            raise KeyError(
                f"The assembly ID '{assembly_id}' is not found"
            )

        # Parse transformations from assembly
        transformations = []
        for transform in selected_assembly:
            matrix = np.array(transform["matrix"]).reshape(4, 4).copy(order = 'C') # order needs to be 'c' otherwise Blender doesn't like it
            chain_ids = np.array(self._file["chainNameList"], dtype="U4")
            affected_chain_ids = chain_ids[transform["chainIndexList"]]
            transformations.append((
                affected_chain_ids.tolist(),
                matrix[:3, :3].tolist(),
                matrix[:3, 3].tolist()
            ))
        
        return transformations
    
    def get_assemblies(self):
        assembly_dict = {}
        for assembly_id in self.list_assemblies():
            assembly_dict[assembly_id] = self.get_transformations(assembly_id)
        
        return assembly_dict