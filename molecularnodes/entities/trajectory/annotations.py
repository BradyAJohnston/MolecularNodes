import typing
from MDAnalysis.core.groups import AtomGroup
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
    selection: str | AtomGroup = "name CA"
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
        if isinstance(params.selection, str):
            # check if selection phrase is valid
            # mda throws exception if invalid
            self.atom_group = universe.select_atoms(params.selection)
        elif isinstance(params.selection, AtomGroup):
            self.atom_group = params.selection
        else:
            raise ValueError(f"Need str or AtomGroup. Got {type(params.selection)}")
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
    selection: typing.Union[str | AtomGroup] = "protein"
    text: str = "COM"

    def defaults(self) -> None:
        params = self.interface
        # show a default pointer
        params.pointer_length = 2
        params.arrow_size = 10

    def validate(self) -> bool:
        params = self.interface
        universe = self.trajectory.universe
        if isinstance(params.selection, str):
            # check if selection phrase is valid
            # mda throws exception if invalid
            self.atom_group = universe.select_atoms(params.selection)
        elif isinstance(params.selection, AtomGroup):
            self.atom_group = params.selection
        else:
            raise ValueError(f"Need str or AtomGroup. Got {type(params.selection)}")
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
    selection1: str | AtomGroup
    selection2: str | AtomGroup
    text1: str = "COM1"
    text2: str = "COM2"

    def validate(self) -> bool:
        params = self.interface
        universe = self.trajectory.universe
        # check selection 1
        if isinstance(params.selection1, str):
            # check if selection phrase is valid
            # mda throws exception if invalid
            self.atom_group1 = universe.select_atoms(params.selection1)
        elif isinstance(params.selection1, AtomGroup):
            self.atom_group1 = params.selection1
        else:
            raise ValueError(f"Need str or AtomGroup. Got {type(params.selection1)}")
        # check selection 2
        if isinstance(params.selection2, str):
            # check if selection phrase is valid
            # mda throws exception if invalid
            self.atom_group2 = universe.select_atoms(params.selection2)
        elif isinstance(params.selection2, AtomGroup):
            self.atom_group2 = params.selection2
        else:
            raise ValueError(f"Need str or AtomGroup. Got {type(params.selection2)}")
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


class CanonicalDihedrals(TrajectoryAnnotation):
    """
    Canonical Dihedrals of a given residue

    Attributes
    ----------
    resid: int
        The residue number

    show_atom_names: bool
        Whether or not to show the individual atom names in the residue

    show_direction: bool
        Whether or not to show the arc indicating the angle direction

    """

    annotation_type = "canonical_dihedrals"

    resid: int
    show_atom_names: bool = True
    show_direction: bool = True

    def defaults(self) -> None:
        params = self.interface
        params.arrow_size = 8

    def validate(self) -> bool:
        params = self.interface
        universe = self.trajectory.universe
        # verify residue in universe - raises exception if invalid
        self.residue = universe.residues[params.resid - 1]
        return True

    def draw(self) -> None:
        params = self.interface
        r = self.residue
        # *_selection() returns None if nothing found
        selections = [
            r.phi_selection(),
            r.psi_selection(),
            r.omega_selection(),
            r.chi1_selection(),
        ]
        symbols = ["ϕ", "ψ", "ω", "χ1"]
        # iterate over the canonical selections
        for si in range(4):
            if selections[si] is None:
                continue
            # iterate over the four atoms in the selection
            for ai in range(4):
                a1 = selections[si].atoms[ai]
                if params.show_atom_names:
                    self.draw_text_3d(a1.position, a1.name)
                if ai == 3:
                    continue
                a2 = selections[si].atoms[ai + 1]
                v1 = a1.position
                v2 = a2.position
                text = None
                if ai == 1:
                    text = f"{symbols[si]} = {selections[si].dihedral.value():1.2f} °"
                    if params.show_direction:
                        center = (v1 + v2) / 2
                        radius = 0.1 * self.distance(v1, v2)  # 10% of distance
                        # draw arc indicating direction of the angle
                        self.draw_circle_3d(
                            center, radius, (v2 - v1), angle=270, c_arrow=True
                        )
                self.draw_line_3d(v1, v2, mid_text=text)
