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
    selection: str = "name CA"
    show_resid: bool = False
    show_segid: bool = False

    def defaults(self) -> None:
        params = self.interface
        # set the default text size and text color
        params.text_size = 16
        params.text_color = (1, 1, 1, 1)

    def validate(self) -> bool:
        params = self.interface
        universe = self.trajectory.universe
        # check if selection phrase is valid - mda throws exception if invalid
        self.atom_group = universe.select_atoms(params.selection)
        return True

    def draw(self) -> None:
        params = self.interface
        # iterate over each atom in the atom group
        for atom in self.atom_group:
            text = atom.name
            if params.show_resid:
                text += f"|res {atom.resid}, {atom.resname}"
            if params.show_segid:
                text += f"|seg {atom.segid}"
            self.draw_text_3d(atom.position, text)


class COM(TrajectoryAnnotation):
    """
    Center-of-Mass Trajectory Annotation

    This annotation shows the center-of-mass of a selection.
    This annotation can be added using the 'add_com' method of the
    trajectory's annotation manager:
        trajectory.annotations.add_com(...)

    Attributes
    ----------
    selection: str
        MDAnalysis selection phrase to select the atom group
    text: str
        Text do display at the center-of-mass

    """

    annotation_type = "com"  # required annotation type

    # annotation inputs
    selection: str = "protein"
    text: str = "COM"

    def defaults(self) -> None:
        params = self.interface
        # show a default pointer
        params.pointer_length = 2
        params.arrow_size = 10

    def validate(self) -> bool:
        universe = self.trajectory.universe
        # check if selection phrase is valid - mda throws exception if invalid
        self.atom_group = universe.select_atoms(self.interface.selection)
        return True

    def draw(self) -> None:
        com = self.atom_group.center_of_mass()
        self.draw_text_3d(com, self.interface.text)


class COMDistance(TrajectoryAnnotation):
    """
    Distance between Center-of-Masses Trajectory Annotation

    This annotation shows the distance between center-of-masses of two selections.
    This annotation can be added using the 'add_com_distance' method of the
    trajectory's annotation manager:
        trajectory.annotations.add_com_distance(...)

    Attributes
    ----------
    selection1: str
        MDAnalysis selection phrase to select the first atom group
    selection2: str
        MDAnalysis selection phrase to select the second atom group
    text1: str
        Text do display at the first center-of-mass
    text2: str
        Text do display at the second center-of-mass

    """

    annotation_type = "com_distance"  # required annotation type

    # annotation inputs
    selection1: str
    selection2: str
    text1: str = "COM1"
    text2: str = "COM2"

    def validate(self) -> bool:
        params = self.interface
        universe = self.trajectory.universe
        # check if selection phrases are valid - exception thrown if invalid
        self.atom_group1 = universe.select_atoms(params.selection1)
        self.atom_group2 = universe.select_atoms(params.selection2)
        return True

    def draw(self) -> None:
        params = self.interface
        # calculate the two center-of-masses
        com1 = self.atom_group1.center_of_mass()
        com2 = self.atom_group2.center_of_mass()
        # calculate distance between coms
        distance = self.distance(com1, com2)
        # draw line with arrow ends showing distance in between
        self.draw_line_3d(
            v1=com1,
            v2=com2,
            v1_text=params.text1,
            v2_text=params.text2,
            mid_text=f"{distance:1.2f} Å",
            v1_arrow=True,
            v2_arrow=True,
        )
