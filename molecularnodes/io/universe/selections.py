import bpy
import MDAnalysis as mda
import numpy.typing as npt
import numpy as np
from bpy.props import StringProperty, BoolProperty
from .handlers import _update_universes


class Selection:
    def __init__(
        self, universe: mda.Universe, selection_str, name, updating=True, periodic=True
    ):
        self.selection_str: str = selection_str
        self.periodic: bool = periodic
        self.updating: bool = updating
        self.universe: mda.Universe = universe
        self.message: str = ""
        self.name: str = name
        self.cleanup: bool = True
        self.ag = universe.select_atoms(
            selection_str, updating=updating, periodic=periodic
        )
        self.mask_array = self._ag_to_mask()

    def _ag_to_mask(self) -> npt.NDArray[np.bool_]:
        "Return a 1D boolean mask for the Universe atoms that are in the Selection's AtomGroup."
        return np.isin(self.universe.atoms.ix, self.ag.ix).astype(bool)

    def change_selection(
        self,
        selection_str: str,
        name: str,
        updating: bool = True,
        periodic: bool = True,
    ) -> None:
        "Change the current AtomGroup, using the parent universe and creating a new selection with the given `selectrion_str`"
        self.name = name
        self.periodic = periodic
        self.updating = updating
        self.selection_str = selection_str
        try:
            self.ag = self.universe.select_atoms(
                selection_str, updating=updating, periodic=periodic
            )
            self.message = ""
        except Exception as e:
            self.message = str(e)
            print(e)

    def to_mask(self) -> npt.NDArray[np.bool_]:
        "Returns the selection as a 1D numpy boolean mask. If updating=True, recomputes selection."
        if self.updating:
            self.mask_array = self._ag_to_mask()
        return self.mask_array


class TrajectorySelectionItem(bpy.types.PropertyGroup):
    """Group of properties for custom selections for MDAnalysis import."""

    name: StringProperty(  # type: ignore
        name="Name",
        description="Name of the attribute on the mesh",
        default="custom_selection",
        update=_update_universes,
    )

    selection_str: StringProperty(  # type: ignore
        name="Selection",
        description="Selection to be applied, written in the MDAnalysis selection language",
        default="name CA",
        update=_update_universes,
    )

    updating: BoolProperty(  # type: ignore
        name="Updating",
        description="Recalculate the selection on scene frame change",
        default=True,
        update=_update_universes,
    )

    periodic: BoolProperty(  # type: ignore
        name="Periodic",
        description="For geometric selections, whether to account for atoms in different periodic images when searching",
        default=True,
        update=_update_universes,
    )
    message: StringProperty(  # type: ignore
        name="Message",
        description="Message to report back from `universe.select_atoms()`",
        default="",
    )


class MN_UL_TrajectorySelectionListUI(bpy.types.UIList):
    """UI List"""

    def draw_item(
        self, context, layout, data, item, icon, active_data, active_propname, index
    ):
        custom_icon = "VIS_SEL_11"

        if self.layout_type in {"DEFAULT", "COMPACT"}:
            row = layout.row()
            if item.message != "":
                custom_icon = "ERROR"
                row.alert = True

            row.prop(item, "name", text="", emboss=False)
            row.prop(item, "updating", icon_only=True, icon="FILE_REFRESH")
            row.prop(item, "periodic", icon_only=True, icon="CUBE")

        elif self.layout_type in {"GRID"}:
            layout.alignment = "CENTER"
            layout.label(text="", icon=custom_icon)


class MN_OT_Universe_Selection_Add(bpy.types.Operator):
    "Add a new custom selection to a trajectory"

    bl_idname = "mn.universe_selection_add"
    bl_label = "+"
    bl_description = "Add a new boolean attribute for the given MDA selection string"

    def execute(self, context):
        bob = context.active_object
        bob.mn_universe_selections.add()
        i = int(len(bob.mn_universe_selections) - 1)
        bob.mn_universe_selections[i].name = f"selection_{i + 1}"
        bob.mn["list_index"] = i
        _update_universes(self, context)

        return {"FINISHED"}


class MN_OT_Universe_Selection_Delete(bpy.types.Operator):
    bl_idname = "mda.delete_item"
    bl_label = "-"
    bl_description = "Delete the given boolean selection from the universe"

    @classmethod
    def poll(cls, context):
        return context.active_object.mn_universe_selections

    def execute(self, context):
        bob = context.active_object
        index = bob.mn.universe_selection_index

        sel_list = bob.mn_universe_selections
        sel_list.remove(index)
        bob.mn.universe_selection_index = len(sel_list) - 1
        _update_universes(self, context)

        return {"FINISHED"}


CLASSSES = [
    TrajectorySelectionItem,  # has to be registered before the others to work properly
    MN_UL_TrajectorySelectionListUI,
    MN_OT_Universe_Selection_Add,
    MN_OT_Universe_Selection_Delete,
]
