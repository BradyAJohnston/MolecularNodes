"""
Generate Gaussian cube files of aspirin (PubChem CID 2244) with PySCF.

RHF/6-31G on the PubChem 3D conformer in `tests/data/2244.sdf`. Writes the HOMO,
LUMO, electron density and molecular electrostatic potential.

    uv run --no-project --with pyscf python docs/dev/pyscf_cubes.py OUT_DIR [N]

`tests/data/2244_homo.cube.gz` is the HOMO at N=24, gzipped.
"""

import sys
from pathlib import Path
from pyscf import gto, scf
from pyscf.tools import cubegen

out = Path(sys.argv[1])
n = int(sys.argv[2]) if len(sys.argv) > 2 else 80
sdf = Path(__file__).parents[2] / "tests/data/2244.sdf"

lines = sdf.read_text().splitlines()
n_atoms = int(lines[3][:3])
atoms = [
    f"{line[31:34].strip()} {line[0:10]} {line[10:20]} {line[20:30]}"
    for line in lines[4 : 4 + n_atoms]
]
mol = gto.M(atom="; ".join(atoms), basis="6-31g", unit="Angstrom")
mf = scf.RHF(mol).run()
homo = mol.nelectron // 2 - 1
dm = mf.make_rdm1()

out.mkdir(parents=True, exist_ok=True)
kwargs = dict(nx=n, ny=n, nz=n, margin=4.0)
cubegen.orbital(mol, str(out / "2244_homo.cube"), mf.mo_coeff[:, homo], **kwargs)
cubegen.orbital(mol, str(out / "2244_lumo.cube"), mf.mo_coeff[:, homo + 1], **kwargs)
cubegen.density(mol, str(out / "2244_density.cube"), dm, **kwargs)
cubegen.mep(mol, str(out / "2244_mep.cube"), dm, **kwargs)
