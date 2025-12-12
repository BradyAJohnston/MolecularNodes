import io
from collections import deque
from uuid import uuid1
import bpy
import matplotlib
import matplotlib.pyplot as plt
from MDAnalysis.core.groups import AtomGroup
from PIL import Image
import molecularnodes as mn

matplotlib.use("Agg")

# Annotations are auto-registered
# Unregister previous class if any while debugging / iterating code
if hasattr(
    mn.entities.trajectory.TrajectoryAnnotationManager,
    "add_streaming_com_distance_chart",
):
    mn.entities.trajectory.TrajectoryAnnotationManager.unregister_type(
        "streaming_com_distance_chart"
    )


class StreamingComDistanceChart(mn.entities.trajectory.TrajectoryAnnotation):
    annotation_type = "streaming_com_distance_chart"

    selection1: str | AtomGroup = "resid 1"
    selection2: str | AtomGroup = "resid 129"
    text1: str = ""
    text2: str = ""

    max_values: int = 100
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
        params = self.interface
        params.line_arrow_size = 0.07
        params.text_offset_y = 5
        self._steps = deque(maxlen=params.max_values)
        self._distances = deque(maxlen=params.max_values)
        self._prev_step = None
        self._chart_name = str(uuid1())

    def validate(self, input_name=None):
        params = self.interface
        u = self.trajectory.universe

        if not isinstance(self.trajectory, mn.entities.trajectory.StreamingTrajectory):
            raise ValueError("This annotation requires a Streaming Trajectory")

        label1 = params.text1
        label2 = params.text2
        # check selection 1
        if isinstance(params.selection1, str):
            # check if selection phrase is valid
            # mda throws exception if invalid
            self._ag1 = u.select_atoms(params.selection1)
            if not label1:
                label1 = params.selection1
        elif isinstance(params.selection1, AtomGroup):
            self._ag1 = params.selection1
            if not label1:
                label1 = "selection1"
        else:
            raise ValueError(f"Need str or AtomGroup. Got {type(params.selection1)}")
        # check selection 2
        if isinstance(params.selection2, str):
            # check if selection phrase is valid
            # mda throws exception if invalid
            self._ag2 = u.select_atoms(params.selection2)
            if not label2:
                label2 = params.selection2
        elif isinstance(params.selection2, AtomGroup):
            self._ag2 = params.selection2
            if not label2:
                label2 = "selection2"
        else:
            raise ValueError(f"Need str or AtomGroup. Got {type(params.selection2)}")
        # setup text
        self._label1 = label1
        self._label2 = label2
        self._title = f"Distance between '{label1}' and '{label2}'"
        # reset values
        if input_name in ("selection1", "selection2", "max_values"):
            self._steps = deque(maxlen=params.max_values)
            self._distances = deque(maxlen=params.max_values)
            self._prev_step = None
        elif input_name in ("text1", "text2"):
            self._prev_step = None
        return True

    def draw(self) -> None:
        params = self.interface
        u = self.trajectory.universe
        ts = u.trajectory.ts

        # chart_name = "streaming_chart_1"
        chart_name = self._chart_name
        step = ts.data["step"]
        # calculate distance and add series points
        com1 = self._ag1.center_of_mass()
        com2 = self._ag2.center_of_mass()
        # calculate distance between coms
        distance = self.distance(com1, com2)

        # create plot only on a new step
        if self._prev_step == step:
            chart_image = bpy.data.images[chart_name]
        else:
            # update plot points
            self._distances.append(distance)
            self._steps.append(step)
            # create plot and save to buffer
            buf = io.BytesIO()
            plt.plot(self._steps, self._distances)
            plt.ylabel("Distance (Å)")
            plt.xlabel("Step")
            plt.title(self._title)
            plt.grid(True)
            plt.savefig(buf, format="png", dpi=100)
            plt.close()
            buf.seek(0)
            # use PIL to convert buffer to pixels for blender
            pil_image = Image.open(buf)
            chart_image = self.pil_image_to_bpy_image(pil_image, chart_name)

        # draw line with arrow ends showing distance in between
        self.draw_line_3d(
            v1=com1,
            v2=com2,
            v1_text=self._label1,
            v2_text=self._label2,
            mid_text=f"{distance:1.2f} Å",
            v1_arrow=True,
            v2_arrow=True,
        )
        # draw bpy image
        self.draw_bpy_image(params.location, chart_image, params.scale)
        # update previous time step
        self._prev_step = step
