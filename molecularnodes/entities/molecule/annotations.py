from biotite.structure import AtomArrayStack
from ...annotations.base import BaseAnnotation
from ...annotations.manager import BaseAnnotationManager
from ..annotations import Label2D, Label3D
from ..base import EntityType


class MoleculeAnnotation(BaseAnnotation):
    """
    Base class for a Molecule Annotation

    All molecule annotations should derive from this base class and implement
    the 'draw' method. All derived classes will have access to the molecule
    instance (self.molecule) and all the annotation inputs and common params
    via self.interface.<property>

    An optional 'defaults' method can be provided to set default values
    to the annotation.

    """

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        # Auto-register any sub classes with the annotation manager
        MoleculeAnnotationManager.register(cls)

    def __init__(self, molecule):
        # Allow access to the molecule entity within the annotations
        self.molecule = molecule
        super().__init__()


class MoleculeAnnotationManager(BaseAnnotationManager):
    """
    Annotation Manager for Molecule Entity

    """

    _entity_type = EntityType.MOLECULE
    _classes = {}  # Entity class specific annotation classes

    def __init__(self, entity):
        super().__init__(entity)
        self._interfaces = {}  # Entity instance specific annotation interfaces


class MoleculeInfo(MoleculeAnnotation):
    """
    Molecule Info Annotation

    Attributes
    ----------
    location: tuple[float, float]
        Normalized coordinates (0.0 - 1.0) of the postion in viewport / render

    show_models: bool
        Whether or not to show the number of models in the molecule

    show_atoms: bool
        Whether or not to show the number of atoms in the molecule

    custom_text: str
        Any custom text to add at the end of the annotation

    """

    annotation_type = "molecule_info"

    location: tuple[float, float] = (0.025, 0.05)
    show_models: bool = True
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
        molecule = self.molecule
        text = ""
        if params.show_models:
            text = f"Models: {molecule.n_models}"
        if params.show_atoms:
            array = molecule.array
            if isinstance(array, AtomArrayStack):
                text = text + f"|Atoms[0]: {len(array[0])}"
            else:
                text = text + f"|Atoms: {len(array)}"
        if params.custom_text != "":
            text = text + "|" + params.custom_text
        # Draw text at normalized coordinates wrt viewport / render
        self.draw_text_2d_norm(params.location, text)


class Label2D(MoleculeAnnotation, Label2D):
    """Common Label2D Annotation for all entities"""


class Label3D(MoleculeAnnotation, Label3D):
    """Common Label3D Annotation for all entities"""
