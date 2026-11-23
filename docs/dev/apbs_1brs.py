"""
Regenerate `tests/data/1brs_AD_apbs.dx.gz`: the APBS electrostatic potential of one
barnase-barstar pair from PDB 1BRS (barnase chain A, barstar chain D).

    APBS=/path/to/apbs uv run docs/dev/apbs_1brs.py [--dime 65 65 65] [--out FILE]

Steps: download 1BRS from RCSB, keep the protein atoms of chains A and D (chain D is
the barstar in contact with barnase A), assign AMBER charges and radii with pdb2pqr
(PROPKA protonation at pH 7, run through `uvx`), then solve the linearised
Poisson-Boltzmann equation with APBS using the input below. Values are in kT/e.

Made with pdb2pqr 3.7.1 and APBS 3.4.1 (the Linux release binary from
https://github.com/Electrostatics/apbs/releases). Neither is a Molecular Nodes
dependency. Point `APBS` at the binary; for the release zip also set
`LD_LIBRARY_PATH` to its `lib` folder.

The coarse and fine grid lengths are the ones pdb2pqr's `--apbs-input` chose for this
pair (the molecule's extent plus 20 Å for the fine grid). The fine grid is 65 points
per axis (about 1 Å spacing) instead of pdb2pqr's 161 x 129 x 129, and values are
rounded to 3 significant digits, which keeps the gzipped file under the 500 kB
limit of the large-file pre-commit hook.
"""

import argparse
import gzip
import os
import subprocess
import tempfile
import urllib.request
from pathlib import Path

APBS_INPUT = """read
    mol pqr 1brs_AD.pqr
end
elec
    mg-auto
    dime {dime}
    cglen 89.0936 76.8077 77.7818
    fglen 72.4080 65.1810 65.7540
    cgcent mol 1
    fgcent mol 1
    mol 1
    lpbe
    bcfl sdh
    pdie 2.0000
    sdie 78.5400
    srfm smol
    chgm spl2
    sdens 10.00
    srad 1.40
    swin 0.30
    temp 298.15
    calcenergy total
    calcforce no
    write pot dx 1brs_AD_apbs
end
print elecEnergy 1 end
quit
"""

REPO = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--dime", nargs=3, type=int, default=[65, 65, 65])
    parser.add_argument(
        "--out", type=Path, default=REPO / "tests/data/1brs_AD_apbs.dx.gz"
    )
    parser.add_argument("--digits", type=int, default=3)
    parser.add_argument("--apbs", default=os.environ.get("APBS", "apbs"))
    args = parser.parse_args()

    work = Path(tempfile.mkdtemp(prefix="apbs_1brs_"))
    with urllib.request.urlopen("https://files.rcsb.org/download/1BRS.pdb") as r:
        lines = r.read().decode().splitlines()
    # one barnase-barstar pair, protein atoms only (waters are HETATM)
    pair = [line for line in lines if line.startswith("CRYST1")]
    for chain in "AD":
        pair += [
            line for line in lines if line.startswith("ATOM") and line[21] == chain
        ]
        pair.append("TER")
    (work / "1brs_AD.pdb").write_text("\n".join(pair + ["END"]) + "\n")

    subprocess.run(
        [
            "uvx",
            "--from",
            "pdb2pqr==3.7.1",
            "pdb2pqr30",
            "--ff=AMBER",
            "--titration-state-method=propka",
            "--with-ph=7.0",
            "--drop-water",
            "1brs_AD.pdb",
            "1brs_AD.pqr",
        ],
        cwd=work,
        check=True,
    )
    dime = " ".join(str(d) for d in args.dime)
    (work / "1brs_AD.in").write_text(APBS_INPUT.format(dime=dime))
    subprocess.run([args.apbs, "1brs_AD.in"], cwd=work, check=True)

    # round the values (APBS writes 7 significant digits, which barely compress)
    lines, in_data = [], False
    for line in (work / "1brs_AD_apbs.dx").read_text().splitlines():
        if line.startswith(("object", "attribute", "component", "#")):
            in_data = "data follows" in line
            lines.append(line)
        elif in_data and line.strip():
            lines.append(" ".join(f"{float(v):.{args.digits}g}" for v in line.split()))
        else:
            lines.append(line)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.out, "wt", compresslevel=9) as dst:
        dst.write("\n".join(lines) + "\n")
    print(f"wrote {args.out} ({args.out.stat().st_size / 1e6:.2f} MB); work dir {work}")


if __name__ == "__main__":
    main()
