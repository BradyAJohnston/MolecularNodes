"""
Trajectory selection management for MolecularNodes.

This module manages atom selections for an mda.Universe, coordinating between MDAnalysis
``AtomGroup`` objects, Blender UI properties, and geometry node attributes.

Notes
-----
Architecture
~~~~~~~~~~~~
The ``SelectionManager`` maintains selections through three synchronized layers:

1. **UI Layer**: Blender ``CollectionProperty`` (``mn_trajectory_selections``)
   - Source of truth for selection state
   - User-editable properties (string, updating, periodic)
   - Displays error messages from failed selections

2. **Python Layer**: Cached ``AtomGroup`` objects in ``SelectionManager.atomgroups``
   - Created from UI items via :meth:`SelectionManager.ui_item_to_ag`
   - Recreated when selection strings change
   - Removed when UI items are deleted

3. **Geometry Layer**: Boolean attributes on Blender mesh
   - Stored via :meth:`SelectionManager.ag_to_attribute`
   - Updated for changed or dynamic (``UpdatingAtomGroup``) selections
   - Used by geometry nodes for visual styling

Data Flow
---------
User creates selection → UI item created → :meth:`SelectionManager.update_attributes`
→ ``AtomGroup`` created → Boolean attribute stored → Geometry nodes use attribute

Classes
-------
SelectionError
    Exception for selection operation failures.
FrozenUpdates
    Context manager to freeze updates during bulk UI changes.
SelectionManager
    Main manager coordinating all selection operations.

See Also
--------
molecularnodes.entities.trajectory.base.Trajectory : Uses SelectionManager.
molecularnodes.ui.props.TrajectorySelectionItem : UI property definition.
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
    """Context manager to temporarily freeze SelectionManager updates.

    Prevents :meth:`update_attributes` from running during multiple UI changes.

    Parameters
    ----------
    manager : SelectionManager
        The manager to freeze.

    See Also
    --------
    SelectionManager.from_string : Uses this context manager internally.
    SelectionManager.from_atomgroup : Uses this context manager internally.
    """

    def __init__(self, manager: SelectionManager):
        self.manager = manager

    def __enter__(self):
        self.manager._is_frozen = True

    def __exit__(self, type, value, traceback):
        self.manager._is_frozen = False


class SelectionManager:
    """Manages atom selections for a trajectory.

    Coordinates between MDAnalysis ``AtomGroup`` objects, Blender UI properties, and
    geometry node attributes. Selections are stored as boolean attributes
    on the trajectory object for use in geometry nodes.

    The ``CollectionProperty`` is the 'source of truth' for managing selections for the
    trajectory. If an ``AtomGroup`` doesn't have a matching UI Item in the collection
    property, it will be discarded. New ``AtomGroup`` objects are created for new UI Items.

    The collection is registered and available under ``mn_trajectory_selections`` on an
    object inside of Blender. It can be accessed on this class via :attr:`ui_items` and
    individual items via ``self.ui_items.get('name')``.

    ```python
    bpy.types.Object.mn_trajectory_selections = CollectionProperty(
        type=props.TrajectorySelectionItem
    )
    ```

    Parameters
    ----------
    trajectory : Trajectory
        Parent trajectory object.

    Attributes
    ----------
    atomgroups : dict[str, mda.AtomGroup]
        Cached ``AtomGroup`` objects keyed by selection name.
    ui_index : IntObjectMNProperty
        Property descriptor for current UI selection index.
    """

    ui_index = IntObjectMNProperty("trajectory_selection_index")

    def __init__(self, trajectory: Trajectory):
        self.trajectory = trajectory
        self.atomgroups: dict[str, mda.AtomGroup] = {}
        self._is_frozen: bool = False

    @property
    def ui_items(self) -> bpy.types.bpy_prop_collection_idprop:
        """Blender UI property collection that represent the selections for the Trajectory.

        Returns
        -------
        bpy.types.bpy_prop_collection_idprop
            Collection of :class:`TrajectorySelectionItem` objects stored on the Blender object
            as ``mn_trajectory_selections``.
        """
        return self.object.mn_trajectory_selections  # type: ignore

    @property
    def object(self) -> bpy.types.Object:
        """Blender object associated with the Trajectory and SelectionManager.

        Returns
        -------
        bpy.types.Object
            The Blender object representing the trajectory mesh.
        """
        return self.trajectory.object

    @property
    def universe(self) -> mda.Universe:
        """MDAnalysis Universe for the Trajectory and SelectionManager.

        Returns
        -------
        mda.Universe
            The MDAnalysis Universe containing topology and trajectory data.
        """
        return self.trajectory.universe

    def ag_to_attribute(self, ag: mda.AtomGroup, name: str) -> None:
        """Convert and store an AtomGroup as a boolean attribute.

        Converts an ``AtomGroup`` to a boolean mask into the original ``Universe`` that would
        return the selected atoms in the ``AtomGroup``. This array is then stored as a boolean
        attribute on the mesh that represents the ``Universe`` inside of Blender.

        Parameters
        ----------
        ag : mda.AtomGroup
            The atom group to convert.
        name : str
            Name for the attribute.

        See Also
        --------
        _ag_to_bool : Helper function that performs the AtomGroup to boolean conversion.
        update_attributes : Calls this method to sync selections to geometry attributes.
        """
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
        """Create a selection from an MDAnalysis selection string.

        This uses the MDAnalysis selection language to create an ``AtomGroup`` and stores the
        selection of which atoms are in the ``AtomGroup`` as a boolean attribute on the mesh
        inside of Blender.

        Parameters
        ----------
        string : str
            MDAnalysis selection string (e.g., ``"protein"``, ``"resid 1-10"``).
        updating : bool, default=True
            If ``True``, selection potentially updates each frame if required (e.g., distance-based
            selections). If ``False``, creates a static selection.
        periodic : bool, default=True
            Consider periodic boundary conditions for geometric selections (e.g., ``"around"``).
        name : str, optional
            Name for the selection, used as the attribute name when storing on the mesh.
            Auto-generated if not provided via :meth:`_unique_selection_name`.

        Returns
        -------
        TrajectorySelectionItem
            The created UI item for the selection.

        See Also
        --------
        from_atomgroup : Create selection from pre-existing AtomGroup.
        update_attributes : Called after item creation to generate the AtomGroup.
        _unique_selection_name : Generates unique names when not provided.
        """
        with FrozenUpdates(self):
            item: TrajectorySelectionItem = self.ui_items.add()
            item.name = name if name else self._unique_selection_name()
            item.string = string
            item.updating = updating
            item.periodic = periodic

        self.update_attributes()

        return item

    def _unique_selection_name(self) -> str:
        """Generate a unique selection name.

        Attempts to not clash with existing attribute and selection names.

        Returns
        -------
        str
            Unique name like ``"selection"``, ``"selection_1"``, etc.

        See Also
        --------
        _unique_aname : Utility function that generates unique attribute names.
        """
        return _unique_aname(self.object, "selection")

    def ag_is_updating(self, atomgroup: mda.AtomGroup) -> bool:
        """Check if an AtomGroup is an UpdatingAtomGroup.

        ``UpdatingAtomGroup`` objects recalculate their members each frame based on
        geometric criteria (e.g., distance-based selections).

        Parameters
        ----------
        atomgroup : mda.AtomGroup
            The atom group to check.

        Returns
        -------
        bool
            ``True`` if the ``AtomGroup`` updates dynamically, ``False`` if static.

        Notes
        -----
        Uses class name comparison since ``UpdatingAtomGroup`` is a subclass of ``AtomGroup``.
        """
        return atomgroup.__class__.__name__ == "UpdatingAtomGroup"

    def from_atomgroup(
        self,
        atomgroup: mda.AtomGroup,
        *,
        name: str | None = None,
    ) -> TrajectorySelectionItem:
        """Create a selection from an existing MDAnalysis AtomGroup.

        Create a selection on the Trajectory from an already created ``AtomGroup`` rather than
        just using a string selection input. The selection string displayed is non-editable
        in the GUI.

        Parameters
        ----------
        atomgroup : mda.AtomGroup
            Pre-existing ``AtomGroup`` (static or updating).
        name : str, optional
            Name for the selection. Auto-generated if not provided via :meth:`_unique_selection_name`.

        Returns
        -------
        TrajectorySelectionItem
            The created UI item for the selection with ``item.from_atomgroup = True``.

        Notes
        -----
        Sets ``item.from_atomgroup = True`` to prevent string editing in UI.
        The string representation is stored for display purposes only.

        See Also
        --------
        from_string : Create selection from MDAnalysis selection string.
        ag_to_attribute : Called to immediately store the selection as an attribute.
        """
        with FrozenUpdates(self):
            item: TrajectorySelectionItem = self.ui_items.add()  # type: ignore
            item.name = name if name else self._unique_selection_name()
            self.atomgroups[item.name] = atomgroup
            item.from_atomgroup = True
            ag_as_string = str(atomgroup)
            item.string = ag_as_string
            item.previous_string = ag_as_string

        self.ag_to_attribute(atomgroup, item.name)
        return item

    def ui_item_to_ag(self, item: TrajectorySelectionItem) -> mda.AtomGroup:
        """Generate an ``mda.AtomGroup`` from a ``TrajectorySelectionItem``.

        Uses the item's ``string``, ``updating``, and ``periodic`` properties to create
        the corresponding ``AtomGroup`` from the trajectory's ``Universe``.

        Parameters
        ----------
        item : TrajectorySelectionItem
            The UI item containing selection parameters (``string``, ``updating``, ``periodic``).

        Returns
        -------
        mda.AtomGroup
            ``AtomGroup`` (or ``UpdatingAtomGroup``) created from the item's parameters.

        See Also
        --------
        update_attributes : Calls this method to create missing AtomGroups.
        """
        return self.universe.select_atoms(
            item.string, updating=item.updating, periodic=item.periodic
        )

    def update_attributes(self) -> None:
        """Synchronize UI items, AtomGroups, and named attributes.

        This is the core update method called when selections change. The following steps
        are carried out:

        1. Creates missing ``AtomGroup`` objects for UI items
        2. Removes orphaned ``AtomGroup`` objects with no matching UI item
        3. Create new ``AtomGroup`` objects when selection strings change on the UI item
        4. Update the named attributes on the mesh for updating or new selections

        Any errors in creation are stored as ``item.message`` which will be reflected in the
        UI with a warning and the error message.

        Notes
        -----
        Skipped when manager is frozen via :class:`FrozenUpdates` context when creating new UI items.

        See Also
        --------
        ui_item_to_ag : Creates AtomGroups from UI items.
        ag_to_attribute : Stores AtomGroups as geometry attributes.
        from_string : Calls this after creating UI item.
        """
        if self._is_frozen:
            return

        # first we iterate through UI items to ensure they have a corresponding atomgroup
        # stored for creating attributes from
        for item in self.ui_items:
            if item.name not in self.atomgroups:
                try:
                    self.atomgroups[item.name] = self.ui_item_to_ag(item)
                except Exception as e:
                    item.message = str(e)
                    continue

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

            if item.updating and item.from_atomgroup:
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
        """Remove a selection by name or index.

        Cleans up the UI item, cached ``AtomGroup``, and geometry attribute. Silently
        handles cases where attribute or ``AtomGroup`` don't exist.

        Parameters
        ----------
        value : int or str
            Selection name (str) or index (int) in ``ui_items`` collection.

        Raises
        ------
        ValueError
            If name not found in ``ui_items`` or value is neither int nor str.

        See Also
        --------
        update_attributes : Automatically removes orphaned AtomGroups.
        """
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

    def get(self, name: str) -> TrajectorySelectionItem | None:
        """Try and get a selection UI Item by name.

        Parameters
        ----------
        name : str
            Name of the UI item to retrieve.

        Returns
        -------
        TrajectorySelectionItem or None
            The matching UI item, or ``None`` if no match was found.
        """
        return self.ui_items.get(name)

    def __len__(self) -> int:
        """Return the number of selections.

        Enables use of ``len(manager)`` to get selection count.

        Returns
        -------
        int
            Number of UI items representing the number of selections.

        See Also
        --------
        ui_items : The underlying collection being measured.
        """
        return len(self.ui_items)

    def __getitem__(self, name: str) -> Any:
        """Get a UI selection item by name.

        Enables dictionary-style access: ``manager['selection_name']``.

        Parameters
        ----------
        name : str
            The selection name.

        Returns
        -------
        TrajectorySelectionItem
            The UI item for the selection.

        Raises
        ------
        KeyError
            If no selection with the given name exists.

        See Also
        --------
        get : Returns None instead of raising KeyError if name not found.
        ui_items : The underlying collection being accessed.
        """
        return self.ui_items[name]
