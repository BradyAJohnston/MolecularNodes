import io
from uuid import uuid1
import addon_utils
import bpy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis import helix_analysis as hel
from PIL import Image

addon_id = "bl_ext.blender_org.molecularnodes"
is_enabled, is_loaded = addon_utils.check(addon_id)
if is_enabled and is_loaded:
    import bl_ext.blender_org.molecularnodes as mn  # type: ignore
else:
    import molecularnodes as mn

matplotlib.use("Agg")

# Annotations are auto-registered
# Unregister previous class if any while debugging / iterating code
if hasattr(
    mn.entities.trajectory.TrajectoryAnnotationManager,
    "add_helix_analysis",
):
    mn.entities.trajectory.TrajectoryAnnotationManager.unregister_type("helix_analysis")


class TrajectoryHelixAnalysis(mn.entities.trajectory.TrajectoryAnnotation):
    annotation_type = "helix_analysis"

    # selection has to be atleast 9 residues
    selection: str = "name CA and resnum 161-187"
    chart: list[str] = [
        "local_twists",
        "local_bends",
        "local_heights",
        "local_nres_per_turn",
        "local_origins",
        "local_axes",
        "local_helix_directions",
        "local_screw_angles",
        "global_axis",
    ]

    show_local_axes: bool = False
    show_global_axis: bool = False
    local_axes_length: float = 1.0
    global_axis_length: float = -15.0

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

        # run analysis initially or when selection input changes
        if input_name in (None, "selection") or self._h is None:
            # check selection
            if isinstance(params.selection, str):
                # check if selection phrase is valid
                # mda throws exception if invalid
                self._ag = u.select_atoms(params.selection)
            else:
                self._h = None
                raise ValueError(f"Need str. Got {type(params.selection)}")

            # save current frame
            current_frame = u.trajectory.frame
            # Helix analysis using HELNAL
            # From: https://userguide.mdanalysis.org/stable/examples/analysis/structure/helanal.html
            self._h = hel.HELANAL(u, select=params.selection).run()
            # restore current frame
            u.trajectory[current_frame]
            if not self._h.results.summary:
                self._h = None
                raise ValueError("Invalid selection phrase for helix analysis")
            self._prev_frame = None
        # recreate chart if chart type changes
        if input_name in (None, "chart"):
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
            results = getattr(self._h.results, params.chart)
            mean_values = results.mean(axis=1)
            plt.plot(mean_values)
            if isinstance(mean_values[frame], np.ndarray):
                for v in mean_values[frame]:
                    plt.scatter(frame, v, c="r", s=100, zorder=1)
            else:
                plt.scatter(frame, mean_values[frame], c="r", s=100, zorder=1)
            plt.title(params.chart)
            plt.ylabel("Average value")
            plt.xlabel("Frame")
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
        # show local axis if enabled
        if params.show_local_axes:
            local_origins = self._h.results.local_origins[frame]
            local_axes = self._h.results.local_axes[frame]
            for i, local_axis in enumerate(local_axes):
                self.draw_line_3d(
                    v1=local_origins[i],
                    v2=(local_origins[i] + (local_axis * params.local_axes_length)),
                    v2_arrow=True,
                )
        # show global axis if enabled
        if params.show_global_axis:
            local_origins = self._h.results.local_origins[frame]
            global_axis = self._h.results.global_axis[frame]
            self.draw_line_3d(
                v1=local_origins[0],
                v2=(local_origins[0] + (global_axis * params.global_axis_length)),
                v2_arrow=True,
                overrides={"line_arrow_size": 0.1},
            )
        # update previous frame
        self._prev_frame = frame
