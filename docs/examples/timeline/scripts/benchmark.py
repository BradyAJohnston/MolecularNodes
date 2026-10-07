"""Reproducible end-to-end timings for a turntable of replicated ubiquitin.

After ProteinMotion's ``examples/benchmark.py``: ``copies`` translated copies
of 1UBQ, one full turn, timed. Not a claim about peak performance.
"""

import argparse
import json
import platform
import time
from _common import QUICK, make_canvas, output_dir, spin
import molecularnodes as mn


def benchmark(copies: int, style: str, seconds: float, fps: int):
    canvas = make_canvas(fps=fps)
    edge = int(round(copies ** (1 / 3))) or 1
    molecules = []
    for i in range(copies):
        mol = mn.Molecule.fetch("1UBQ")
        mol.add_style(style)
        mol.object.location = (
            3.8 * (i % edge),
            3.8 * ((i // edge) % edge),
            3.8 * (i // (edge * edge)),
        )
        molecules.append(mol)
    canvas.look_at(
        molecules[0].object
        if copies == 1
        else sum((m.get_view() for m in molecules[1:]), molecules[0].get_view()),
        viewpoint="default",
        margin=0.1,
    )
    with canvas.timeline() as t:
        t.play(*(spin(m, 360) for m in molecules), run_time=seconds, easing="linear")
    started = time.perf_counter()
    t.render(output_dir() / f"benchmark_{copies}_{style}.mp4")
    elapsed = time.perf_counter() - started
    frames = t.frame_end - t.start + 1
    return {
        "copies": copies,
        "style": style,
        "frames": frames,
        "seconds": round(elapsed, 2),
        "frames_per_second": round(frames / elapsed, 2),
        "engine": type(canvas.engine).__name__,
        "resolution": list(canvas.resolution),
        "platform": platform.platform(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("out_dir", nargs="?")
    parser.add_argument("--copies", type=int, default=1 if QUICK else 8)
    parser.add_argument("--style", default="cartoon")
    parser.add_argument("--seconds", type=float, default=1.0 if QUICK else 4.0)
    args = parser.parse_args()
    report = benchmark(args.copies, args.style, args.seconds, fps=24)
    (output_dir() / "benchmark.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
