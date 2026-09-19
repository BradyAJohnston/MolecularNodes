# Node-group asset "Charge" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from .fallback_float import FallbackFloat


class Charge(AssetGeometryGroup):
    """
    Read the `charge` attribute from the geometry. The attribute is only present when the file or topology provided charges (partial charges from a simulation topology, formal charges from a structure file's charge column or `pdbx_formal_charge`). Outputs 0.0 when the geometry has no `charge` attribute.

    Parameters
    ----------
    index : InputInteger
        Index

    Inputs
    ------
    i.index : IntegerSocket
        Index

    Outputs
    -------
    o.charge : FloatSocket
        Read the `charge` attribute from the geometry, only present when the file or topology provided charges. 0.0 if the geometry has no `charge` attribute
    """

    _name = "Charge"
    _asset_name = "Charge"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Read the `charge` attribute from the geometry. The attribute is only present when the file or topology provided charges (partial charges from a simulation topology, formal charges from a structure file's charge column or `pdbx_formal_charge`). Outputs 0.0 when the geometry has no `charge` attribute.",
        "node_tool_idname": "geometry.charge",
    }

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        charge: FloatSocket
        """Read the `charge` attribute from the geometry, only present when the file or topology provided charges. 0.0 if the geometry has no `charge` attribute"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        index: InputInteger = 0,
    ):
        super().__init__(**{"Index": index})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        index = tree.inputs.integer("Index", 0, min_value=0, default_input="INDEX")
        charge = tree.outputs.float(
            "charge",
            description="Read the `charge` attribute from the geometry, only present when the file or topology provided charges. 0.0 if the geometry has no `charge` attribute",
        )

        FallbackFloat(name="charge").o.value.point.at(index) >> charge


ASSET = Charge

ASSET_METADATA = {
    "description": "Read the `charge` attribute from the geometry. The attribute is only present when the file or topology provided charges (partial charges from a simulation topology, formal charges from a structure file's charge column or `pdbx_formal_charge`). Outputs 0.0 when the geometry has no `charge` attribute.",
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
