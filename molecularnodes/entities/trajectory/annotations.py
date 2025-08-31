import os
import typing
from mathutils import Vector
from MDAnalysis.core.groups import AtomGroup
from ...annotations.base import BaseAnnotation
from ...annotations.manager import BaseAnnotationManager
from ..annotations import Label2D, Label3D
from ..base import EntityType


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

    _entity_type = EntityType.MD
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

    def defaults(self) -> None:
        params = self.interface
        params.arrow_size = 10

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

    def _get_start_dv(self, dihedral: list) -> Vector:
        """Get the direction vector along which to start the arc"""
        v0 = Vector(dihedral.atoms[0].position)
        v1 = Vector(dihedral.atoms[1].position)
        v2 = Vector(dihedral.atoms[2].position)
        # find two vectors on the first plane (between 0, 1, 2)
        v10 = v0 - v1
        v20 = v0 - v2
        # the normal of the plane on which the arc is drawn
        normal = v2 - v1
        # (v10 x normal) x v20 is a vector in the first plane that is
        # perpendicular to the normal - get the unit vector
        dv = v10.cross(normal).cross(v20)
        dv.normalize()
        return dv

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
                    angle = selections[si].dihedral.value()
                    text = f"{symbols[si]} = {angle:1.2f} °"
                    if params.show_direction:
                        center = (v1 + v2) / 2
                        radius = 0.1 * self.distance(v1, v2)  # 10% of distance
                        # draw arc indicating direction of the angle
                        start_dv = self._get_start_dv(selections[si])
                        self.draw_circle_3d(
                            center,
                            radius,
                            (v2 - v1),
                            angle=angle,
                            start_dv=start_dv,
                            c_arrow=True,
                        )
                self.draw_line_3d(v1, v2, mid_text=text)


class UniverseInfo(TrajectoryAnnotation):
    """
    Universe Info Trajectory Annotation

    Attributes
    ----------
    location: tuple[float, float]
        Normalized coordinates (0.0 - 1.0) of the postion in viewport / render

    show_frame: bool
        Whether or not to show the frame number

    show_topology: bool
        Whether or not to show the topology filename

    show_trajectory: bool
        Whether or not to show the trajectory filename

    show_atoms: bool
        Whether or not to show the number of atoms

    custom_text: str
        Any custom text to add at the end of the annotation

    """

    annotation_type = "universe_info"

    location: tuple[float, float] = (0.025, 0.05)
    show_frame: bool = True
    show_topology: bool = True
    show_trajectory: bool = True
    show_atoms: bool = True
    custom_text: str = ""

    def defaults(self) -> None:
        params = self.interface
        params.text_align = "left"

    def validate(self) -> bool:
        params = self.interface
        x, y = params.location
        if (not 0 <= x <= 1) or (not 0 <= y <= 1):
            raise ValueError("Normalized coordinates should lie between 0 and 1")
        return True

    def draw(self) -> None:
        params = self.interface
        u = self.trajectory.universe
        text = ""
        if params.show_frame:
            text = f"Frame : {u.trajectory.frame} / {u.trajectory.n_frames - 1}"
        if params.show_topology and u.filename:
            text = text + "|Topology : " + os.path.basename(u.filename)
        if params.show_trajectory and u.trajectory.filename:
            text = text + "|Trajectory : " + os.path.basename(u.trajectory.filename)
        if params.show_atoms:
            text = text + "|Atoms : " + str(u.trajectory.n_atoms)
        if params.custom_text != "":
            text = text + "|" + params.custom_text
        # Draw text at normalized coordinates wrt viewport / render
        self.draw_text_2d_norm(params.location, text)


class Label2D(TrajectoryAnnotation, Label2D):
    """Common Label2D Annotation for all entities"""


class Label3D(TrajectoryAnnotation, Label3D):
    """Common Label3D Annotation for all entities"""
