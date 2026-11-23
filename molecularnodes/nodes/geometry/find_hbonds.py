# Node-group asset "Find HBonds" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry
from ._shared.mn_units import MNUnits
from .atomic_number import AtomicNumber


class FindHBonds(AssetGeometryGroup):
    """
    Find hydrogen bonds between the atoms, as fields on the hydrogens

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry to find hydrogen bonds in
    selection : InputBoolean
        Selection of atoms considered as hydrogens, donors and acceptors
    distance : InputFloat
        Maximum distance between the donor and acceptor in Angstroms
    angle : InputFloat
        Minimum donor-hydrogen-acceptor angle

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry to find hydrogen bonds in
    i.selection : BooleanSocket
        Selection of atoms considered as hydrogens, donors and acceptors
    i.distance : FloatSocket
        Maximum distance between the donor and acceptor in Angstroms
    i.angle : FloatSocket
        Minimum donor-hydrogen-acceptor angle

    Outputs
    -------
    o.is_hbonded : BooleanSocket
        True for each hydrogen in a hydrogen bond
    o.acceptor_index : IntegerSocket
        Index of the acceptor a hydrogen bonds to, -1 for other atoms
    """

    _name = "Find HBonds"
    _asset_name = "Find HBonds"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Find hydrogen bonds between the atoms, as fields on the hydrogens"
    }

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry to find hydrogen bonds in"""
        selection: BooleanSocket
        """Selection of atoms considered as hydrogens, donors and acceptors"""
        distance: FloatSocket
        """Maximum distance between the donor and acceptor in Angstroms"""
        angle: FloatSocket
        """Minimum donor-hydrogen-acceptor angle"""

    class _Outputs(SocketAccessor):
        is_hbonded: BooleanSocket
        """True for each hydrogen in a hydrogen bond"""
        acceptor_index: IntegerSocket
        """Index of the acceptor a hydrogen bonds to, -1 for other atoms"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        distance: InputFloat = 3.0,
        angle: InputFloat = 5 * math.pi / 6,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Distance": distance,
                "Angle": angle,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry to find hydrogen bonds in"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms considered as hydrogens, donors and acceptors",
            hide_value=True,
        )
        distance = tree.inputs.float(
            "Distance",
            3.0,
            description="Maximum distance between the donor and acceptor in Angstroms",
            min_value=0.0,
            max_value=10.0,
        )
        angle = tree.inputs.float(
            "Angle",
            5 * math.pi / 6,
            description="Minimum donor-hydrogen-acceptor angle",
            min_value=0.0,
            max_value=math.pi,
            subtype="ANGLE",
        )
        is_hbonded = tree.outputs.boolean(
            "Is HBonded", description="True for each hydrogen in a hydrogen bond"
        )
        acceptor_index = tree.outputs.integer(
            "Acceptor Index",
            description="Index of the acceptor a hydrogen bonds to, -1 for other atoms",
        )

        atomic_number = AtomicNumber()
        boolean_math = g.Compare.integer.equal(
            atomic_number, 7
        ).o.result | g.Compare.integer.equal(atomic_number, 8)
        capture = g.CaptureAttribute.point(geometry=atoms)
        atom_index = capture.items.integer(g.Index(), "Atom Index")
        position = g.Position()
        with g.Frame("Donors"):
            separate_geometry = g.SeparateGeometry.point(
                capture.o.geometry, atomic_number > 1
            )
            sample_index = g.SampleIndex.point.integer(
                separate_geometry.o.selection,
                atom_index.output,
                g.SampleNearest.point(separate_geometry.o.selection),
            )
            evaluate_at_index = position.o.position.point.at(sample_index)
            _string = g.String(
                "Each hydrogen is paired with the nearest heavy atom within 1.2 A, as MDAnalysis does, rather than by bonds, as hydrogens often have none (N-terminal hydrogens, deuterium). Deuterium has no atomic number so reads as 0, and counts as hydrogen."
            )
            boolean_math_1 = (
                selection
                & (atomic_number <= 1)
                & (
                    position.o.position.distance(evaluate_at_index)
                    <= MNUnits(value=1.2).o.angstrom
                )
            )
            boolean_math_2 = (
                boolean_math_1
                & boolean_math.point.at(sample_index)
                & selection.point.at(sample_index)
            )
        with g.Frame("Acceptors"):
            mn_units = MNUnits(value=distance)
            math_1 = angle.cos()
            separate_geometry_1 = g.SeparateGeometry.point(
                capture.o.geometry, boolean_math & selection
            )
            _string_1 = g.String(
                "The nearest N or O to a probe point past the hydrogen, along donor to hydrogen, so the donor itself is never the nearest. Probes at 1.0, 0.5 and 1.5 A are tried in turn, and the first acceptor within the distance and angle cutoffs is kept. Each hydrogen bonds to at most one acceptor."
            )
            vector_math = (position.o.position - evaluate_at_index).normalize()
            vector_math_1 = (evaluate_at_index - position).normalize()
            sample_nearest = g.SampleNearest.point(
                separate_geometry_1.o.selection,
                position.o.position + vector_math * MNUnits(value=1.5).o.angstrom,
            )
            sample_index_1 = g.SampleIndex.point.integer(
                separate_geometry_1.o.selection, atom_index.output, sample_nearest
            )
            evaluate_at_index_1 = position.o.position.point.at(sample_index_1)
            boolean_math_3 = (
                boolean_math_2
                & g.Compare.integer.not_equal(sample_index_1, sample_index)
                & (
                    evaluate_at_index.distance(evaluate_at_index_1)
                    <= mn_units.o.angstrom
                )
            )
            boolean_math_4 = boolean_math_3 & (
                vector_math_1.dot((evaluate_at_index_1 - position).normalize())
                <= math_1
            )
            sample_nearest_1 = g.SampleNearest.point(
                separate_geometry_1.o.selection,
                position.o.position + vector_math * MNUnits(value=0.5).o.angstrom,
            )
            sample_index_2 = g.SampleIndex.point.integer(
                separate_geometry_1.o.selection, atom_index.output, sample_nearest_1
            )
            evaluate_at_index_2 = position.o.position.point.at(sample_index_2)
            boolean_math_5 = (
                boolean_math_2
                & g.Compare.integer.not_equal(sample_index_2, sample_index)
                & (
                    evaluate_at_index.distance(evaluate_at_index_2)
                    <= mn_units.o.angstrom
                )
            )
            boolean_math_6 = boolean_math_5 & (
                vector_math_1.dot((evaluate_at_index_2 - position).normalize())
                <= math_1
            )
            switch = boolean_math_6.switch.integer(
                boolean_math_4.switch.integer(g.Integer(-1), sample_index_1),
                sample_index_2,
            )
            sample_nearest_2 = g.SampleNearest.point(
                separate_geometry_1.o.selection,
                position.o.position + vector_math * MNUnits(value=1.0).o.angstrom,
            )
            sample_index_3 = g.SampleIndex.point.integer(
                separate_geometry_1.o.selection, atom_index.output, sample_nearest_2
            )
            evaluate_at_index_3 = position.o.position.point.at(sample_index_3)
            boolean_math_7 = (
                boolean_math_2
                & g.Compare.integer.not_equal(sample_index_3, sample_index)
                & (
                    evaluate_at_index.distance(evaluate_at_index_3)
                    <= mn_units.o.angstrom
                )
            )
            boolean_math_8 = boolean_math_7 & (
                vector_math_1.dot((evaluate_at_index_3 - position).normalize())
                <= math_1
            )
            switch_1 = boolean_math_8.switch.integer(switch, sample_index_3)
        (switch_1 >= 0) >> is_hbonded

        switch_1 >> acceptor_index


ASSET = FindHBonds

ASSET_METADATA = {
    "description": "Find hydrogen bonds between the atoms, as fields on the hydrogens",
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
