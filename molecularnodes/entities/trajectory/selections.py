"""
Trajectory selection management for MolecularNodes.

This module provides classes for managing atom selections in molecular dynamics
trajectories. It bridges MDAnalysis AtomGroups with Blender's geometry node system,
allowing dynamic selections that update during trajectory playback.

Classes
-------
Selection
    Wraps an MDAnalysis AtomGroup and syncs it with Blender UI properties
    and geometry node attributes.
SelectionManager
    Manages all selections for a trajectory, coordinating between Python objects,
    Blender UI, and geometry nodes.

Notes
-----
Selections can be either static (computed once) or updating (recalculated every frame).
The selection state is stored as boolean mask arrays that are exposed as geometry
node attributes for use in Blender's node-based rendering system.
"""

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
    """
    Represents an atom selection for a trajectory in Blender.

    A Selection wraps an MDAnalysis AtomGroup and provides synchronization between
    the MDAnalysis selection and Blender UI properties. It manages the selection string,
    updating behavior, and conversion to boolean masks for geometry node attributes.

    Parameters
    ----------
    trajectory : Trajectory
        The parent Trajectory object this selection belongs to.
    name : str, optional
        Unique name for this selection. Default is "selection_0".

    Attributes
    ----------
    trajectory : Trajectory
        Reference to the parent trajectory.
    mask_array : npt.NDArray[np.bool_]
        Cached boolean mask representing which atoms are selected.
    """

    def __init__(self, trajectory: Trajectory, name: str = "selection_0"):
        self._ag: mda.AtomGroup | None = None
        self._uuid: str = str(uuid1())
        self._name = name
        self._current_string: str = ""
        self.trajectory = trajectory

    def add_selection_property(self, string: str = "all", updating=True, periodic=True):
        """
        Initialize the selection with a UI property and create the atom group.

        Creates a Blender UI property item for this selection and sets up the
        MDAnalysis AtomGroup based on the provided selection string.

        Parameters
        ----------
        string : str, optional
            MDAnalysis selection string. Default is "all".
        updating : bool, optional
            Whether the selection should update dynamically during trajectory playback.
            Default is True.
        periodic : bool, optional
            Whether to consider periodic boundary conditions in the selection.
            Default is True.
        """
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
        """
        Get the name of this selection.

        Returns
        -------
        str
            The unique name identifying this selection.
        """
        return self._name

    @property
    def ui_item(self):
        """
        Get the Blender UI property item for this selection.

        Returns
        -------
        TrajectorySelectionItem
            The Blender property group containing UI-accessible properties
            for this selection (string, updating, periodic, message, etc.).
        """
        return self.trajectory.selections.items[self.name]

    @property
    def string(self) -> str:
        """
        Get the MDAnalysis selection string.

        Returns
        -------
        str
            The current selection string (e.g., "name CA", "resid 1:10").
        """
        return self.ui_item.string

    @string.setter
    def string(self, string: str) -> None:
        """
        Set the MDAnalysis selection string and update the atom group.

        Updates the selection string in the UI property and attempts to
        create a new AtomGroup. If successful, updates the selection mask
        in Blender. If an error occurs, stores the error message.

        Parameters
        ----------
        string : str
            MDAnalysis selection string to apply.
        """
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
        """
        Get whether periodic boundary conditions are considered.

        Returns
        -------
        bool
            True if periodic boundary conditions are enabled for this selection.
        """
        return self.ui_item.periodic

    @periodic.setter
    def periodic(self, periodic: bool) -> None:
        """
        Set whether to consider periodic boundary conditions.

        Parameters
        ----------
        periodic : bool
            If True, periodic boundary conditions will be considered in the selection.
        """
        if self.ui_item.periodic == periodic:
            return
        self.ui_item.periodic = periodic

    @property
    def updating(self) -> bool:
        """
        Get whether this selection updates dynamically during playback.

        Returns
        -------
        bool
            True if the selection recalculates every frame during trajectory playback.
        """
        return self.ui_item.updating

    @updating.setter
    def updating(self, updating: bool) -> None:
        """
        Set whether this selection should update dynamically during playback.

        Parameters
        ----------
        updating : bool
            If True, the selection will recalculate every frame. If False,
            the selection is static and uses the mask from the initial frame.
        """
        if self.ui_item.updating == updating:
            return
        self.ui_item.updating = updating

    @property
    def message(self) -> str:
        """
        Get the current status or error message for this selection.

        Returns
        -------
        str
            Error message if selection parsing failed, empty string otherwise.
        """
        return self.ui_item.message

    @message.setter
    def message(self, message: str) -> None:
        """
        Set the status or error message for this selection.

        Parameters
        ----------
        message : str
            Message to display to the user, typically an error from selection parsing.
        """
        self.ui_item.message = message

    @property
    def immutable(self) -> bool:
        """
        Get whether this selection can be modified by the user.

        Returns
        -------
        bool
            True if the selection is locked from user modification.
        """
        return self.ui_item.immutable

    @immutable.setter
    def immutable(self, immutable: bool) -> None:
        """
        Set whether this selection can be modified by the user.

        Parameters
        ----------
        immutable : bool
            If True, prevents the user from modifying this selection in the UI.
        """
        self.ui_item.immutable = immutable

    def set_atom_group(self, string: str) -> None:
        """
        Create or update the MDAnalysis AtomGroup from a selection string.

        Parses the selection string using MDAnalysis and creates an AtomGroup
        (static or updating depending on the updating property). Updates the
        cached mask array if successful.

        Parameters
        ----------
        string : str
            MDAnalysis selection string to parse.

        Notes
        -----
        If the selection string is invalid, stores the error message and prints
        a warning but does not raise an exception.
        """
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
        """
        Convert the AtomGroup to a boolean mask array.

        Creates a boolean mask where True indicates atoms in the selection and
        False indicates atoms not in the selection. The mask is aligned with the
        trajectory's Universe atom ordering.

        Returns
        -------
        npt.NDArray[np.bool_]
            1D boolean array with length equal to the number of atoms in the
            Universe, where True indicates the atom is in this selection.
        """
        return np.isin(self.trajectory.universe.atoms.ix, self._ag.ix).astype(bool)

    def set_selection(self) -> None:
        """
        Store the selection mask as a named attribute in Blender.

        Converts the selection to a boolean mask and stores it as a geometry
        node attribute on the trajectory object. Only applies if updating=True.

        Notes
        -----
        This method is called after the selection string is updated to sync
        the selection state with the Blender geometry nodes.
        """
        if not self.updating:
            return
        self.trajectory.store_named_attribute(
            self.to_mask(), name=self.name, atype=db.AttributeTypes.BOOLEAN
        )

    def to_mask(self) -> npt.NDArray[np.bool_]:
        """
        Get the selection as a boolean mask array.

        Returns the cached mask array. If updating=True, recomputes the mask
        from the current frame's AtomGroup state before returning.

        Returns
        -------
        npt.NDArray[np.bool_]
            1D boolean array indicating which atoms are selected.
        """
        if self.updating:
            self.mask_array = self._ag_to_mask()
        return self.mask_array

    @classmethod
    def from_atomgroup(cls, trajectory, atomgroup: mda.AtomGroup, name: str = ""):
        """
        Create a Selection from an existing MDAnalysis AtomGroup.

        This factory method creates a Selection object from an existing AtomGroup,
        automatically detecting whether it's an UpdatingAtomGroup and extracting
        the selection string and periodic settings if available.

        Parameters
        ----------
        trajectory : Trajectory
            The parent trajectory this selection belongs to.
        atomgroup : mda.AtomGroup
            An MDAnalysis AtomGroup or UpdatingAtomGroup to wrap.
        name : str, optional
            Name for the selection. If empty, uses the selection string or
            a default name based on atom count.

        Returns
        -------
        Selection
            A new Selection object wrapping the provided AtomGroup.

        Notes
        -----
        - If atomgroup is an UpdatingAtomGroup, extracts the selection string
          and periodic settings from it.
        - The selection is marked as immutable to prevent user modification.
        - The selection is automatically registered with the trajectory's
          SelectionManager.
        """
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
    """
    Manages all atom selections for a trajectory.

    The SelectionManager coordinates between Selection objects (Python side),
    Blender UI properties, and geometry node attributes. It provides methods
    to create, retrieve, and remove selections, and maintains synchronization
    between the different representation layers.

    Parameters
    ----------
    trajectory : Trajectory
        The parent trajectory object this manager belongs to.

    Attributes
    ----------
    trajectory : Trajectory
        Reference to the parent trajectory.
    """

    def __init__(self, trajectory: "Trajectory"):
        self.trajectory = trajectory
        self._selections: Dict[str, Selection] = {}

    def add(
        self,
        string: str,
        name: str = "selection_0",
        updating: bool = True,
        periodic: bool = False,
    ) -> Selection:
        """
        Create and register a new selection from a selection string.

        Creates a new Selection object with the specified parameters and
        registers it with both the manager's internal dictionary and the
        Blender UI property collection.

        Parameters
        ----------
        string : str
            MDAnalysis selection string (e.g., "name CA", "resid 1:10").
        name : str, optional
            Unique name for this selection. Default is "selection_0".
        updating : bool, optional
            Whether the selection updates every frame. Default is True.
        periodic : bool, optional
            Whether to consider periodic boundary conditions. Default is False.

        Returns
        -------
        Selection
            The newly created and registered Selection object.
        """
        sel = self._selections[name] = Selection(self.trajectory, name=name)
        sel.add_selection_property(string=string, updating=updating, periodic=periodic)
        return sel

    def append(self, selection: Selection) -> None:
        """
        Register an existing Selection object with this manager.

        Adds a Selection object to the internal dictionary without creating
        a new UI property (assumes the Selection was created via a classmethod
        that already handled UI property creation).

        Parameters
        ----------
        selection : Selection
            The Selection object to register.
        """
        self._selections[selection.name] = selection

    def from_ui_item(self, item: TrajectorySelectionItem) -> Selection:
        """
        Create a Selection from a Blender UI property item.

        Creates a new Selection based on the properties stored in a Blender
        UI property item, typically used when loading from a saved .blend file.

        Parameters
        ----------
        item : TrajectorySelectionItem
            Blender property group containing selection parameters.

        Returns
        -------
        Selection
            The newly created Selection object.
        """
        return self.add(
            item.string, item.name, updating=item.updating, periodic=item.periodic
        )

    def from_atomgroup(
        self, atomgroup: mda.AtomGroup, name: str = "NewSelection"
    ) -> Selection:
        """
        Create a Selection from an existing MDAnalysis AtomGroup.

        Wraps an AtomGroup in a Selection object and stores the initial mask
        as a Blender geometry attribute. Useful for creating selections from
        programmatically defined AtomGroups.

        Parameters
        ----------
        atomgroup : mda.AtomGroup
            MDAnalysis AtomGroup or UpdatingAtomGroup to wrap.
        name : str, optional
            Name for the new selection. Default is "NewSelection".

        Returns
        -------
        Selection
            The newly created Selection object.

        See Also
        --------
        Selection.from_atomgroup : The underlying classmethod used to create the Selection.
        """
        sel = Selection.from_atomgroup(
            trajectory=self.trajectory, atomgroup=atomgroup, name=name
        )
        self._selections[sel.name] = sel
        self.trajectory.store_named_attribute(
            sel.to_mask(), name=sel.name, atype=db.AttributeTypes.BOOLEAN
        )
        return sel

    def remove(self, name: str) -> None:
        """
        Remove a selection from the manager and clean up all references.

        Removes the selection from the internal dictionary, the Blender UI
        property collection, and the geometry node attribute. Handles cases
        where the selection may be partially registered.

        Parameters
        ----------
        name : str
            Name of the selection to remove.

        Notes
        -----
        Silently ignores KeyError or AttributeError if the selection was not
        fully registered in all locations.
        """
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
        """
        Retrieve a Selection by name.

        Parameters
        ----------
        name : str
            Name of the selection to retrieve.

        Returns
        -------
        Selection
            The Selection object with the given name.

        Raises
        ------
        KeyError
            If no selection with the given name exists.
        """
        return self._selections[name]

    @property
    def items(self) -> bpy.types.CollectionProperty:
        """
        Get the Blender UI property collection for all selections.

        Returns
        -------
        bpy.types.CollectionProperty
            The Blender property collection containing TrajectorySelectionItem
            properties for this trajectory.
        """
        return self.trajectory.object.mn_trajectory_selections

    @property
    def index(self) -> int:
        """
        Get the currently selected index in the UI list.

        Returns
        -------
        int
            The index of the currently selected selection in the UI list.
        """
        return self.trajectory.object.mn["list_index"]

    @index.setter
    def index(self, value: int) -> None:
        """
        Set the currently selected index in the UI list.

        Parameters
        ----------
        value : int
            The index to select in the UI list.
        """
        self.trajectory.object.mn["list_index"] = value

    def __len__(self) -> int:
        """
        Get the number of selections managed by this manager.

        Returns
        -------
        int
            The number of registered Selection objects.
        """
        return len(self._selections)
