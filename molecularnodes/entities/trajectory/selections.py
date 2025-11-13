"""
Trajectory selection management for MolecularNodes.
"""

from __future__ import annotations
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ...ui.props import TrajectorySelectionItem
    from .base import Trajectory
import bpy
import databpy as db
import MDAnalysis as mda
from ..utilities import IntObjectMNProperty, _unique_aname
from .helpers import _ag_to_bool


class SelectionError(Exception):
    """Exception raised when selection operations fail."""

    pass


class FrozenUpdates:
    def __init__(self, manager: SelectionManager):
        self.manager = manager

    def __enter__(self):
        self.manager._is_frozen = True

    def __exit__(self, type, value, traceback):
        self.manager._is_frozen = False


class SelectionManager:
    ui_index = IntObjectMNProperty("trajectory_selection_index")

    def __init__(self, trajectory: Trajectory):
        self.trajectory = trajectory
        self.atomgroups: dict[str, mda.AtomGroup] = {}
        self._is_frozen: bool = False

    @property
    def ui_items(self) -> bpy.types.bpy_prop_collection_idprop:
        return self.object.mn_trajectory_selections  # type: ignore

    @property
    def object(self) -> bpy.types.Object:
        return self.trajectory.object

    @property
    def universe(self) -> mda.Universe:
        return self.trajectory.universe

    def ag_to_attribute(self, ag: mda.AtomGroup, name: str) -> None:
        self.trajectory.store_named_attribute(
            data=_ag_to_bool(ag), name=name, atype=db.AttributeTypes.BOOLEAN
        )

    def from_string(
        self,
        string: str,
        *,
        updating: bool = True,
        periodic: bool = True,
        name: str | None = None,
    ) -> TrajectorySelectionItem:
        with FrozenUpdates(self):
            item: TrajectorySelectionItem = self.ui_items.add()
            item.name = name if name else self._unique_selection_name()
            item.string = string
            item.updating = updating
            item.periodic = periodic

        self.update_attributes()

        return item

    def _unique_selection_name(self) -> str:
        return _unique_aname(self.object, "selection")

    def ag_is_updating(self, atomgroup: mda.AtomGroup) -> bool:
        return atomgroup.__class__.__name__ == "UpdatingAtomGroup"

    def from_atomgroup(
        self,
        atomgroup: mda.AtomGroup,
        *,
        name: str | None = None,
        immutable: bool = False,
    ) -> TrajectorySelectionItem:
        with FrozenUpdates(self):
            item: TrajectorySelectionItem = self.ui_items.add()  # type: ignore
            item.name = name if name else self._unique_selection_name()
            self.atomgroups[item.name] = atomgroup
            item.from_atomgroup = True
            ag_as_string = str(atomgroup)
            item.string = ag_as_string
            item.previous_string = ag_as_string
            item.immutable = immutable

        self.ag_to_attribute(atomgroup, item.name)
        return item

    def ui_item_to_ag(self, item: TrajectorySelectionItem) -> mda.AtomGroup:
        return self.universe.select_atoms(
            item.string, updating=item.updating, periodic=item.periodic
        )

    def update_attributes(self) -> None:
        if self._is_frozen:
            return

        # first we iterate through UI items to ensure they have a corresponding atomgroup
        # stored for creating attributes from
        for item in self.ui_items:
            if item.name not in self.atomgroups:
                self.atomgroups[item.name] = self.ui_item_to_ag(item)

        # then we iterate through the stored AtomGroups. If they don't have a corresponding
        # UI item we remove the atomgroup from the dictionary so we can stop caring about it
        # if an item has been created without a name - initialise a unique attribute name
        for key in list(self.atomgroups.keys()):
            selection_has_changed = False
            ag = self.atomgroups[key]

            item: TrajectorySelectionItem = self.ui_items.get(key)  # type: ignore
            if item is None:
                del self.atomgroups[key]
                continue

            if item.name == "":
                item.name = self._unique_selection_name()

            if item.from_atomgroup:
                if self.ag_is_updating(ag):
                    self.ag_to_attribute(ag, item.name)
                continue

            # if the ui selection item can't be successfully created as an AtomGroup the
            # message will be set to the error and this flagged in the UI. Upon successfully
            # creating the AtomGroup we set the messag to "" as a signal everything is OK
            # and continue on with setting the attribute value
            # We also have to track if we need to update the AtomGroup because the user has
            # changed the input string at all. We just track with a "previous_string" which
            # should only be updated if an AtomGroup has successfully been created with the
            # input selection string
            if item.string != item.previous_string or item.message != "":
                try:
                    ag = self.universe.select_atoms(
                        item.string, updating=item.updating, periodic=item.periodic
                    )
                    self.atomgroups[item.name] = ag
                    item.message = ""
                    item.previous_string = item.string
                    selection_has_changed = True
                except Exception as e:
                    item.message = str(e)
                    continue

            if selection_has_changed or (item.updating and self.ag_is_updating(ag)):
                self.ag_to_attribute(ag, item.name)

    def remove(self, value: int | str):
        if isinstance(value, str):
            idx = self.ui_items.find(value)
            if idx == -1:
                raise ValueError("{} not present in UI Items".format(value))
        elif isinstance(value, int):
            idx = value
        else:
            raise ValueError("`value` must be either string or integer")

        item: TrajectorySelectionItem = self.ui_items[idx]

        # we want to remove the named attribute that might be stored, the atomgroup we are
        # tracking with it and also finally the ui_item that corresponds to it all
        try:
            self.trajectory.remove_named_attribute(item.name)
        except db.NamedAttributeError:
            pass
        try:
            del self.atomgroups[item.name]
        except KeyError:
            pass
        self.ui_items.remove(idx)  # type: ignore

    def __len__(self) -> int:
        return len(self.atomgroups)

    def __getitem__(self, name: str) -> Any:
        return self.ui_items[name]
