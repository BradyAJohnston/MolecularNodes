"""
Trajectory selection management for MolecularNodes.
"""

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ...ui.props import TrajectorySelectionItem
    from .base import Trajectory
from uuid import uuid1
import bpy
import databpy as db
import MDAnalysis as mda
from ..utilities import IntObjectMNProperty, _unique_aname
from .helpers import _ag_to_bool


class SelectionError(Exception):
    """Exception raised when selection operations fail."""

    pass


class SelectionManager:
    ui_index = IntObjectMNProperty("trajectory_selection_index")

    def __init__(self, trajectory: Trajectory):
        self.trajectory = trajectory
        self.atomgroups: dict[str, mda.AtomGroup] = {}

    @property
    def ui_items(self) -> bpy.types.CollectionProperty:
        return self.object.mn_trajectory_selections  # type: ignore

    @property
    def object(self) -> bpy.types.Object:
        return self.trajectory.object

    @property
    def universe(self) -> mda.Universe:
        return self.trajectory.universe

    def ag_to_attribute(self, ag: mda.AtomGroup, attribute_name: str) -> None:
        self.trajectory.store_named_attribute(
            data=_ag_to_bool(ag), name=attribute_name, atype=db.AttributeTypes.BOOLEAN
        )

    def add(
        self,
        string: str,
        updating: bool = True,
        periodic: bool = True,
        att_name: str | None = None,
    ) -> TrajectorySelectionItem:
        return self.from_atomgroup(
            self.trajectory.universe.select_atoms(
                string, updating=updating, periodic=periodic
            ),
            att_name=att_name,
        )

    def _unique_selection_name(self) -> str:
        return _unique_aname(self.object, "selection")

    def from_atomgroup(
        self,
        atomgroup: mda.AtomGroup,
        att_name: str | None = None,
        immutable: bool = False,
    ) -> TrajectorySelectionItem:
        item = self.ui_items.add()  # type: ignore
        item.name = str(uuid1())
        item.attribute_name = att_name if att_name else self._unique_selection_name()
        item.immutable = immutable
        self.atomgroups[item.name] = atomgroup
        self.update_attributes()
        return item

    def from_ui_item(self, item: TrajectorySelectionItem) -> None:
        ag = self.trajectory.universe.select_atoms(
            item.string, updating=item.updating, periodic=item.periodic
        )
        self.atomgroups[item.name] = ag

    def update_attributes(self) -> None:
        # first we iterate through UI items to ensure they have a corresponding atomgroup
        # stored for creating attributes from
        for item in self.ui_items:
            sel = self.atomgroups.get(item.name)
            if not sel:
                self.from_ui_item(item)

        # then we iterate through the stored AtomGroups. If they don't have a corresponding
        # UI item we remove the atomgroup from the dictionary so we can stop caring about it
        for name, ag in self.atomgroups.items():
            sel = self.ui_items.get(name)
            if not sel:
                del self.atomgroups[name]
                continue
            if sel.attribute_name == "":
                sel.attribute_name = self._unique_selection_name()

            if sel.previous_string != sel.string or not sel.message == "":
                try:
                    ag = self.universe.select_atoms(
                        sel.string, updating=sel.updating, periodic=sel.periodic
                    )
                    sel.message = ""
                    sel.previous_string = sel.string
                except Exception as e:
                    sel.message = str(e)
                    continue

            if sel.updating:
                self.atomgroups[sel.name] = ag
                sel.previous_string = sel.string

                self.ag_to_attribute(ag, sel.attribute_name)

    def remove(self, value: int | str):
        if isinstance(value, str):
            idx = [i.name for i in self.ui_items].index(value)
        elif isinstance(value, int):
            idx = value
        else:
            raise ValueError("`value` must be either string or integer")
        sel = self.ui_items[idx]  # type: ignore
        try:
            self.trajectory.remove_named_attribute(sel.attribute_name)
        except db.NamedAttributeError:
            pass
        del self.atomgroups[sel.name]
        self.ui_items.remove(idx)  # type: ignore

    def __len__(self) -> int:
        return len(self.atomgroups)
