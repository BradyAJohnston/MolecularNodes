# Node group 'Vector in Angstroms' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputFloat, InputVector
from ..angstrom_to_world import AngstromToWorld


class VectorInAngstroms(CustomGeometryGroup):
    """
    Vector in Angstroms

    Parameters
    ----------
    vector : InputVector
        Vector
    normalize : InputBoolean
        Normalize the vector before first scaling to angstroms
    angstrom : InputFloat
        Angstrom

    Inputs
    ------
    i.vector : VectorSocket
        Vector
    i.normalize : BooleanSocket
        Normalize the vector before first scaling to angstroms
    i.angstrom : FloatSocket
        Angstrom

    Outputs
    -------
    o.vector : VectorSocket
        Vector that has been scaled by the number of input angstroms, optionally normalizing the vector first
    """

    _name = "Vector in Angstroms"
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.vector_in_angstroms"}

    class _Inputs(SocketAccessor):
        vector: VectorSocket
        """Vector"""
        normalize: BooleanSocket
        """Normalize the vector before first scaling to angstroms"""
        angstrom: FloatSocket
        """Angstrom"""

    class _Outputs(SocketAccessor):
        vector: VectorSocket
        """Vector that has been scaled by the number of input angstroms, optionally normalizing the vector first"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        vector: InputVector = None,
        normalize: InputBoolean = True,
        angstrom: InputFloat = 2.5,
    ):
        super().__init__(
            **{"Vector": vector, "Normalize": normalize, "Angstrom": angstrom}
        )

    def _build_group(self, tree):
        vector = tree.inputs.vector(
            "Vector", (0.0, 0.0, 0.0), min_value=-10_000.0, max_value=10_000.0
        )
        normalize = tree.inputs.boolean(
            "Normalize",
            True,
            description="Normalize the vector before first scaling to angstroms",
        )
        angstrom = tree.inputs.float(
            "Angstrom", 2.5, min_value=-10_000.0, max_value=10_000.0
        )
        vector_1 = tree.outputs.vector(
            "Vector",
            description="Vector that has been scaled by the number of input angstroms, optionally normalizing the vector first",
        )

        (
            normalize.switch.vector(vector, vector.normalize())
            * AngstromToWorld(angstrom=angstrom)
            >> vector_1
        )
