import io
from collections import deque
from uuid import uuid1
import bl_ext.blender_org.molecularnodes as mn  # type: ignore
import bpy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.core.groups import AtomGroup
from PIL import Image

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

    selection: str | AtomGroup = "protein"
    periodic: bool = True
    updating: bool = True
    text: str = ""
    show_circles: bool = True

    max_values: int = 100
    location: tuple[float, float] = (0.025, 0.05)
    scale: float = 0.75

    def defaults(self):
        params = self.interface
        plt.rcParams.update(
            {
                # white with alpha = 30%
                "figure.facecolor": (1.0, 1.0, 1.0, 0.3),
                "axes.facecolor": (1.0, 1.0, 1.0, 0.3),
                "savefig.facecolor": (1.0, 1.0, 1.0, 0.3),
            }
        )
        self._steps = deque(maxlen=params.max_values)
        self._rogs = deque(maxlen=params.max_values)
        self._com = None
        self._prev_step = None
        self._chart_name = str(uuid1())
        params.line_width = 4.0
        params.mesh_thickness = 4.0

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
            self._ag = u.select_atoms(
                params.selection, periodic=params.periodic, updating=params.updating
            )
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
        if input_name in ("selection", "periodic", "updating", "max_values"):
            self._steps = deque(maxlen=params.max_values)
            self._rogs = deque(maxlen=params.max_values)
            self._com = None
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
            self._com = self._ag.center_of_mass()
            # get squared distance from center
            ri_sq = (coordinates - self._com) ** 2
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
            self._rog = np.sqrt(rog_sq)
            # update plot points
            self._rogs.append(self._rog)
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

        if params.show_circles and self._com is not None:
            # draw circles with corresponding matplotlib series default colors
            # #1f77b4 (Blue), #ff7f0e (Orange), #2ca02c (Green), #d62728 (Red)
            # overall
            self.draw_circle_3d(
                self._com,
                self._rog[0],
                overrides={
                    "line_color": (0.122, 0.467, 0.706, 1),
                    "mesh_color": (0.122, 0.467, 0.706, 1),
                },
            )
            # along x-axis
            self.draw_circle_3d(
                self._com,
                self._rog[1],
                (0, 1, 0),
                overrides={
                    "line_color": (1.0, 0.498, 0.055, 1),
                    "mesh_color": (1.0, 0.498, 0.055, 1),
                },
            )
            # along y-axis
            self.draw_circle_3d(
                self._com,
                self._rog[2],
                (0, 0, 1),
                overrides={
                    "line_color": (0.173, 0.627, 0.173, 1),
                    "mesh_color": (0.173, 0.627, 0.173, 1),
                },
            )
            # along z-axis
            self.draw_circle_3d(
                self._com,
                self._rog[3],
                (1, 0, 0),
                overrides={
                    "line_color": (0.84, 0.15, 0.16, 1),
                    "mesh_color": (0.84, 0.15, 0.16, 1),
                },
            )

        # draw bpy image
        self.draw_bpy_image(params.location, chart_image, params.scale)
        # update previous time step
        self._prev_step = step
