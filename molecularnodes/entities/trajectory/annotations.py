from ...annotations.base import BaseAnnotation
from ...annotations.manager import BaseAnnotationManager


class TrajectoryAnnotation(BaseAnnotation):
    """
    Base class for a Trajectory Annotation

    All trajectory annotations should derive from this base class and implement
    the 'draw' method. All derived classes will have access to the trajectory
    instance (self.trajectory) and all the annotation inputs and common params
    via self.interface.<property>

    An optional 'defaults' method can be provided to set default values
    to the annotation.

    """

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        # Auto-register any sub classes with the annotation manager
        TrajectoryAnnotationManager.register(cls)

    def __init__(self, trajectory):
        # Allow access to the trajectory entity within the annotations
        self.trajectory = trajectory
        super().__init__()


class TrajectoryAnnotationManager(BaseAnnotationManager):
    """
    Annotation Manager for Trajectory Entity

    """

    _classes = {}  # Entity class specific annotation classes

    def __init__(self, entity):
        super().__init__(entity)
        self._interfaces = {}  # Entity instance specific annotation interfaces


class AtomInfo(TrajectoryAnnotation):
    """
    Atom Info Trajectory Annotation

    This annotation shows the atom info of a selection.
    This annotation can be added using the 'add_atom_info' method of the
    trajectory's annotation manager:
        trajectory.annotations.add_atom_info(...)

    Attributes
    ----------
    selection: str
        MDAnalysis selection phrase to select the atom group
    show_resid: bool
        Whether or not to show the res ID along with the atom name
    show_segid: bool
        Whether or not to show the seg ID along with the atom name

    """

    annotation_type = "atom_info"  # required annotation type

    # annotation inputs
    selection: str = "protein"
    show_resid: bool = False
    show_segid: bool = False

    def defaults(self) -> None:
        params = self.interface
        # set the default color for this annotation to red
        params.text_color = (1, 0, 0, 1)

    def draw(self) -> None:
        params = self.interface
        print(self.trajectory, params.selection, params.show_resid, params.show_segid)
