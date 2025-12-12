import io
from uuid import uuid1
import bpy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.core.groups import AtomGroup
from PIL import Image
import molecularnodes as mn

matplotlib.use("Agg")

# Annotations are auto-registered
# Unregister previous class if any while debugging / iterating code
if hasattr(
    mn.entities.trajectory.TrajectoryAnnotationManager,
    "add_streaming_radii_of_gyration_chart",
):
    mn.entities.trajectory.TrajectoryAnnotationManager.unregister_type(
        "streaming_radii_of_gyration_chart"
    )


class StreamingRadiiOfGyrationChart(mn.entities.trajectory.TrajectoryAnnotation):
    annotation_type = "streaming_radii_of_gyration_chart"

    selection: str | AtomGroup = "resid 1"
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
        self._steps = []
        self._rogs = []
        self._prev_step = None
        self._chart_name = str(uuid1())

    def validate(self, input_name=None):
        params = self.interface
        u = self.trajectory.universe

        if not isinstance(self.trajectory, mn.entities.trajectory.StreamingTrajectory):
            raise ValueError("This annotation requires a Streaming Trajectory")

        label = params.text
        # check selection
        if isinstance(params.selection, str):
            # check if selection phrase is valid
            # mda throws exception if invalid
            self._ag = u.select_atoms(params.selection)
            if not label:
                label = params.selection
        elif isinstance(params.selection, AtomGroup):
            self._ag = params.selection
            if not label:
                label = "selection"
        else:
            raise ValueError(f"Need str or AtomGroup. Got {type(params.selection)}")
        # setup masses
        self._masses = self._ag.masses
        self._total_mass = np.sum(self._masses)
        # setup text
        self._title = f"Radii of Gyration of '{label}'"
        # reset values
        if input_name == "selection":
            self._steps = []
            self._rogs = []
            self._prev_step = None
        elif input_name == "text":
            self._prev_step = None
        return True

    def draw(self) -> None:
        params = self.interface
        u = self.trajectory.universe
        step = u.trajectory.ts.data["step"]

        # create plot only on a new step
        if self._prev_step == step:
            chart_image = bpy.data.images[self._chart_name]
        else:
            # calculate radii of gyration and add series points
            # From: https://userguide.mdanalysis.org/stable/examples/analysis/custom_trajectory_analysis.html#Radius-of-gyration
            coordinates = self._ag.positions
            center_of_mass = self._ag.center_of_mass()
            # get squared distance from center
            ri_sq = (coordinates - center_of_mass) ** 2
            # sum the unweighted positions
            sq = np.sum(ri_sq, axis=1)
            sq_x = np.sum(ri_sq[:, [1, 2]], axis=1)  # sum over y and z
            sq_y = np.sum(ri_sq[:, [0, 2]], axis=1)  # sum over x and z
            sq_z = np.sum(ri_sq[:, [0, 1]], axis=1)  # sum over x and y
            # make into array
            sq_rs = np.array([sq, sq_x, sq_y, sq_z])
            # weight positions
            rog_sq = np.sum(self._masses * sq_rs, axis=1) / self._total_mass
            # square root
            rog = np.sqrt(rog_sq)
            # update plot points
            self._rogs.append(rog)
            self._steps.append(step)
            # create plot and save to buffer
            data = np.array(self._rogs)
            labels = ["all", "x-axis", "y-axis", "z-axis"]
            for i, label in enumerate(labels):
                plt.plot(self._steps, data[:, i], label=label)
            plt.legend(loc="upper left")
            plt.ylabel("Radii (Å)")
            plt.xlabel("Step")
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
        # update previous time step
        self._prev_step = step
