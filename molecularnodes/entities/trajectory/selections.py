from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ...ui.props import TrajectorySelectionItem
    from .base import Trajectory
from typing import Dict
from uuid import uuid1
import bpy
import databpy as db
import MDAnalysis as mda
import numpy as np
import numpy.typing as npt


class Selection:
    def __init__(self, trajectory: Trajectory, name: str = "selection_0"):
        self._ag: mda.AtomGroup | None = None
        self._uuid: str = str(uuid1())
        self._name = name
        self._current_string: str = ""
        self.trajectory = trajectory

    def add_selection_property(self, string: str = "all", updating=True, periodic=True):
        prop = self.trajectory.selections.items.add()
        prop.name = self.name
        prop.uuid = self._uuid
        self.updating = updating
        self.periodic = periodic
        self.set_atom_group(string)
        self.string = string
        self.mask_array = self._ag_to_mask()

    @property
    def name(self) -> str:
        return self._name

    @property
    def ui_item(self):
        return self.trajectory.selections.items[self.name]

    @property
    def string(self) -> str:
        return self.ui_item.string

    @string.setter
    def string(self, string: str) -> None:
        self.ui_item.string = string
        try:
            self.set_atom_group(string)
            self.set_selection()
            self.message = ""
        except Exception as e:
            self.message = str(e)
            print(e)

    @property
    def periodic(self) -> bool:
        return self.ui_item.periodic

    @periodic.setter
    def periodic(self, periodic: bool) -> None:
        if self.ui_item.periodic == periodic:
            return
        self.ui_item.periodic = periodic

    @property
    def updating(self) -> bool:
        return self.ui_item.updating

    @updating.setter
    def updating(self, updating: bool) -> None:
        if self.ui_item.updating == updating:
            return
        self.ui_item.updating = updating

    @property
    def message(self) -> str:
        return self.ui_item.message

    @message.setter
    def message(self, message: str) -> None:
        self.ui_item.message = message

    @property
    def immutable(self) -> bool:
        return self.ui_item.immutable

    @immutable.setter
    def immutable(self, immutable: bool) -> None:
        self.ui_item.immutable = immutable

    def set_atom_group(self, string: str) -> None:
        if self._current_string == string:
            return
        try:
            self._ag = self.trajectory.universe.select_atoms(
                string, updating=self.updating, periodic=self.periodic
            )
            self.mask_array = self._ag_to_mask()
            self._current_string = string
            self.message = ""
        except Exception as e:
            self._current_string = string
            self.message = str(e)
            print(
                str(e)
                + f" in selection: `{self.name}` on object: `{self.trajectory.object.name}`"
            )

    def _ag_to_mask(self) -> npt.NDArray[np.bool_]:
        "Return a 1D boolean mask for the Universe atoms that are in the Selection's AtomGroup."
        return np.isin(self.trajectory.universe.atoms.ix, self._ag.ix).astype(bool)

    def set_selection(self) -> None:
        "Sets the selection in the trajectory"
        if not self.updating:
            return
        self.trajectory.store_named_attribute(
            self.to_mask(), name=self.name, atype=db.AttributeTypes.BOOLEAN
        )

    def to_mask(self) -> npt.NDArray[np.bool_]:
        "Returns the selection as a 1D numpy boolean mask. If updating=True, recomputes selection."
        if self.updating:
            self.mask_array = self._ag_to_mask()
        return self.mask_array

    @classmethod
    def from_atomgroup(cls, trajectory, atomgroup: mda.AtomGroup, name: str = ""):
        "Create a Selection object from an AtomGroup"
        # set default value
        string = f"sel_{atomgroup.n_atoms}_atoms"
        updating = False
        periodic = False

        # if class is an UpdatingAtomGroup
        if atomgroup.__class__.__name__ == "UpdatingAtomGroup":
            updating = True
            # assuming it's a single selection
            # MDA do support `u.select_atoms('index 0', 'around 5 index 0')
            string = atomgroup._selection_strings[0]
            try:
                if atomgroup._selections[0].periodic:
                    periodic = True
            except AttributeError:
                # some selections don't have the periodic attribute
                pass
            except Exception as e:
                print(e)

        if name == "":
            name = string
        sel = cls(trajectory=trajectory, name=name)
        trajectory.selections.append(sel)

        prop = trajectory.object.mn_trajectory_selections.add()
        prop.name = name
        prop.uuid = sel._uuid

        sel._ag = atomgroup
        sel.mask_array = sel._ag_to_mask()
        sel._current_string = name
        sel.updating = updating
        sel.periodic = periodic
        sel.immutable = True
        sel.string = name

        return sel


class SelectionManager:
    def __init__(self, trajectory: "Trajectory"):
        self.trajectory = trajectory
        self._selections: Dict[str, Selection] = {}

    def new(
        self,
        string: str,
        name: str = "selection_0",
        updating: bool = True,
        periodic: bool = False,
    ) -> Selection:
        sel = self._selections[name] = Selection(self.trajectory, name=name)
        sel.add_selection_property(string=string, updating=updating, periodic=periodic)
        return sel

    def append(self, selection: Selection) -> None:
        self._selections[selection.name] = selection

    def from_ui_item(self, item: TrajectorySelectionItem) -> Selection:
        return self.new(
            item.string, item.name, updating=item.updating, periodic=item.periodic
        )

    def from_atomgroup(
        self, atomgroup: mda.AtomGroup, name: str = "NewSelection"
    ) -> Selection:
        sel = Selection.from_atomgroup(
            trajectory=self.trajectory, atomgroup=atomgroup, name=name
        )
        self._selections[sel.name] = sel
        self.trajectory.store_named_attribute(
            sel.to_mask(), name=sel.name, atype=db.AttributeTypes.BOOLEAN
        )
        return sel

    def remove(self, name: str) -> None:
        names = [sel.name for sel in self.items]
        index = names.index(name)
        self.items.remove(index)
        try:
            del self._selections[name]
        except KeyError:
            pass
        try:
            self.trajectory.remove_named_attribute(name)
        except AttributeError:
            pass

    def get(self, name: str) -> Selection:
        return self._selections[name]

    @property
    def items(self) -> bpy.types.CollectionProperty:
        return self.trajectory.object.mn_trajectory_selections

    @property
    def index(self) -> int:
        return self.trajectory.object.mn["list_index"]

    @index.setter
    def index(self, value: int) -> None:
        self.trajectory.object.mn["list_index"] = value

    def __len__(self) -> int:
        return len(self._selections)
