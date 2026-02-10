import io
from uuid import uuid1
import bpy
import matplotlib
import matplotlib.pyplot as plt
from MDAnalysis.analysis import rms
from PIL import Image
import molecularnodes as mn

matplotlib.use("Agg")

# Annotations are auto-registered
# Unregister previous class if any while debugging / iterating code
if hasattr(
    mn.entities.trajectory.TrajectoryAnnotationManager,
    "add_rmsd_chart",
):
    mn.entities.trajectory.TrajectoryAnnotationManager.unregister_type("rmsd_chart")


class TrajectoryRMSDChart(mn.entities.trajectory.TrajectoryAnnotation):
    annotation_type = "rmsd_chart"

    selection: str = "protein"
    periodic: bool = True
    updating: bool = True
    text: str = ""

    location: tuple[float, float] = (0.025, 0.05)
    scale: float = 0.75

    def defaults(self):
        plt.rcParams.update(
            {
                # white with alpha = 30%
                "figure.facecolor": (1.0, 1.0, 1.0, 0.3),
                "axes.facecolor": (1.0, 1.0, 1.0, 0.3),
                "savefig.facecolor": (1.0, 1.0, 1.0, 0.3),
            }
        )
        self._prev_frame = None
        self._chart_name = str(uuid1())

    def validate(self, input_name=None):
        params = self.interface
        u = self.trajectory.universe

        if type(self.trajectory) is not mn.entities.trajectory.Trajectory:
            raise ValueError("This annotation requires a Trajectory entity")

        label = params.text
        # check selection
        if isinstance(params.selection, str):
            # check if selection phrase is valid
            # mda throws exception if invalid
            u.select_atoms(
                params.selection, periodic=params.periodic, updating=params.updating
            )
            if not label:
                label = params.selection
        else:
            raise ValueError(f"Need str. Got {type(params.selection)}")

        # run analysis initially or when selection input changes
        if input_name in (None, "selection", "periodic", "updating"):
            # save current frame
            current_frame = u.trajectory.frame
            # calculate RMSD
            # From: https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/rmsd.html#RMSD-of-a-Universe-with-multiple-selections
            R = rms.RMSD(
                u,  # universe to align
                u,  # reference universe or atomgroup
                select=params.selection,  # group to superimpose and calculate RMSD
                ref_frame=0,
            )  # frame index of the reference
            R.run()
            self._frames = R.results.rmsd[:, 0]
            self._rmsd_values = R.results.rmsd[:, 2]
            # restore current frame
            u.trajectory[current_frame]
        # setup text
        self._title = f"RMSD of '{label}'"
        # reset values
        if input_name in ("selection", "text"):
            self._prev_frame = None
        return True

    def draw(self) -> None:
        params = self.interface
        u = self.trajectory.universe
        frame = u.trajectory.frame

        # create plot only on frame change
        if self._prev_frame == frame:
            chart_image = bpy.data.images[self._chart_name]
        else:
            # create plot and save to buffer
            plt.plot(self._frames, self._rmsd_values)
            plt.scatter(frame, self._rmsd_values[frame], c="r", s=100, zorder=1)
            plt.ylabel("RMSD (Å)")
            plt.xlabel("Frame")
            plt.title(self._title)
            plt.grid(True)
            buf = io.BytesIO()
            plt.savefig(buf, format="png", dpi=100)
            plt.close()
            buf.seek(0)
            # use PIL to convert buffer to pixels for blender
            pil_image = Image.open(buf)
            chart_image = self.pil_image_to_bpy_image(pil_image, self._chart_name)

        # draw bpy image
        self.draw_bpy_image(params.location, chart_image, params.scale)
        # update previous frame
        self._prev_frame = frame
