"""
Trajectory selection management for MolecularNodes.

This module provides classes for managing atom selections in molecular dynamics
trajectories. It bridges MDAnalysis AtomGroups with Blender's geometry node system,
allowing dynamic selections that update during trajectory playback.

Data Flow Architecture
----------------------
PropertyGroups (TrajectorySelectionItem) are the source of truth for persistent data:
    User edits UI → PropertyGroup.update callback → _update_entities() handler
                                                              ↓
                             Selection.sync_from_ui() ← Called by handler
                                                              ↓
                             Selection.set_atom_group() ← Recomputes _ag
                                                              ↓
                             Selection.sync_to_blender_attribute() ← Updates geometry

Classes
-------
SelectionError
    Exception raised when selection operations fail.
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


class SelectionError(Exception):
    """Exception raised when selection operations fail."""

    pass


class AtomGroupMetadata:
    """Metadata extracted from an MDAnalysis AtomGroup."""

    def __init__(self, atomgroup: mda.AtomGroup):
        """
        Initialize metadata from an AtomGroup.

        Parameters
        ----------
        atomgroup : mda.AtomGroup
            The AtomGroup to extract metadata from.
        """
        self._atomgroup = atomgroup
        self._uuid = str(uuid1())

    @property
    def fallback_name(self) -> str:
        return self.string or f"sel_{self._atomgroup.n_atoms}_atoms"

    @property
    def uuid(self) -> str:
        """Get the unique identifier for this selection."""
        return self._uuid

    @property
    def string(self) -> str:
        """Extract the selection string from the AtomGroup."""
        # Assuming it's a single selection
        # MDA does support `u.select_atoms('index 0', 'around 5 index 0')`
        try:
            return self._atomgroup._selection_strings[0]
        except (AttributeError, IndexError):
            return ""

    @property
    def updating(self) -> bool:
        """Determine if the AtomGroup is updating."""
        return self._atomgroup.__class__.__name__ == "UpdatingAtomGroup"

    @property
    def periodic(self) -> bool:
        """Determine if periodic boundary conditions are enabled."""
        if not self.updating:
            return False
        try:
            return self._atomgroup._selections[0].periodic
        except (AttributeError, IndexError):
            # Some selections don't have the periodic attribute
            return False
        except Exception as e:
            print(f"Warning extracting periodic setting: {e}")
            return False


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
        self._ui_item_cache: TrajectorySelectionItem | None = None

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
    def ui_item(self) -> "TrajectorySelectionItem":
        """
        Get the Blender UI property item for this selection.

        Uses caching to avoid repeated lookups. The cache can be invalidated
        by calling invalidate_ui_cache().

        Returns
        -------
        TrajectorySelectionItem
            The Blender property group containing UI-accessible properties
            for this selection (string, updating, periodic, message, etc.).

        Raises
        ------
        RuntimeError
            If the UI property cannot be found, indicating the selection
            was not properly initialized.
        """
        if self._ui_item_cache is None:
            try:
                self._ui_item_cache = self.trajectory.selections.items[self.name]
            except (KeyError, AttributeError) as e:
                raise RuntimeError(
                    f"UI property not found for selection '{self.name}'. "
                    f"Selection may not be properly initialized."
                ) from e
        return self._ui_item_cache

    def invalidate_ui_cache(self) -> None:
        """
        Invalidate the cached UI property reference.

        Call this if the UI item is recreated or the selection is renamed.
        The next access to ui_item will fetch a fresh reference.
        """
        self._ui_item_cache = None

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

        Notes
        -----
        If the selection is marked as immutable, logs a warning and returns
        without making changes.
        """
        if self.immutable:
            print(
                f"Warning: Attempted to modify immutable selection '{self.name}'. "
                "Ignoring change."
            )
            return

        self.ui_item.string = string
        success = self.apply_selection_string(string)
        if success:
            self.sync_to_blender_attribute()

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

    def validate_selection_string(self, string: str) -> tuple[bool, str]:
        """
        Validate a selection string without side effects.

        Tests whether the selection string can be parsed by MDAnalysis
        without actually storing the result.

        Parameters
        ----------
        string : str
            MDAnalysis selection string to validate.

        Returns
        -------
        tuple[bool, str]
            A tuple of (is_valid, error_message). If valid, error_message is empty.

        Examples
        --------
        >>> is_valid, error = sel.validate_selection_string("name CA")
        >>> if is_valid:
        ...     sel.apply_selection_string("name CA")
        """
        try:
            # Just test parsing, don't store
            self.trajectory.universe.select_atoms(
                string, updating=False, periodic=False
            )
            return True, ""
        except Exception as e:
            return False, str(e)

    def apply_selection_string(self, string: str) -> bool:
        """
        Apply a selection string and update the atom group.

        This method validates the selection string and, if valid, creates
        the AtomGroup and updates the cached mask array.

        Parameters
        ----------
        string : str
            MDAnalysis selection string to apply.

        Returns
        -------
        bool
            True if the selection was successfully applied, False otherwise.

        Notes
        -----
        Unlike set_atom_group, this method returns a success status and
        handles errors gracefully by storing them in the message property.
        """
        if self._current_string == string:
            return True

        try:
            self._ag = self.trajectory.universe.select_atoms(
                string, updating=self.updating, periodic=self.periodic
            )
            self.mask_array = self._ag_to_mask()
            self._current_string = string
            self.message = ""
            return True
        except Exception as e:
            error_msg = (
                f"{str(e)} in selection: `{self.name}` "
                f"on object: `{self.trajectory.object.name}`"
            )
            self.message = str(e)
            print(error_msg)
            return False

    def set_atom_group(self, string: str, raise_on_error: bool = False) -> bool:
        """
        Create or update the MDAnalysis AtomGroup from a selection string.

        Parses the selection string using MDAnalysis and creates an AtomGroup
        (static or updating depending on the updating property). Updates the
        cached mask array if successful.

        Parameters
        ----------
        string : str
            MDAnalysis selection string to parse.
        raise_on_error : bool, optional
            If True, raises SelectionError on failure.
            If False, stores error in message and returns False.
            Default is False.

        Returns
        -------
        bool
            True if successful, False if failed (only when raise_on_error=False).

        Raises
        ------
        SelectionError
            If raise_on_error=True and the selection string is invalid.

        Notes
        -----
        This method is maintained for backward compatibility. New code should
        prefer apply_selection_string() for better error handling.
        """
        if self._current_string == string:
            return True

        try:
            self._ag = self.trajectory.universe.select_atoms(
                string, updating=self.updating, periodic=self.periodic
            )
            self.mask_array = self._ag_to_mask()
            self._current_string = string
            self.message = ""
            return True
        except Exception as e:
            error_msg = (
                f"{str(e)} in selection: `{self.name}` "
                f"on object: `{self.trajectory.object.name}`"
            )
            self.message = str(e)
            print(error_msg)

            if raise_on_error:
                raise SelectionError(error_msg) from e
            return False

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

    def sync_to_blender_attribute(self) -> None:
        """
        Synchronize the selection mask to Blender geometry attributes.

        Converts the current selection to a boolean mask and stores it as a
        geometry node attribute on the trajectory object. Only applies if
        updating=True.

        Notes
        -----
        This pushes the Python-side selection state to Blender's geometry system,
        making it available for use in geometry nodes.
        """
        if not self.updating:
            return
        mask = self.to_mask()
        self.trajectory.store_named_attribute(
            mask, name=self.name, atype=db.AttributeTypes.BOOLEAN
        )

    def set_selection(self) -> None:
        """
        Store the selection mask as a named attribute in Blender.

        .. deprecated::
            Use sync_to_blender_attribute() instead for clearer intent.

        Notes
        -----
        This method is maintained for backward compatibility.
        """
        self.sync_to_blender_attribute()

    def sync_from_ui(self) -> None:
        """
        Synchronize selection state from the UI property to computed state.

        Pulls the current state from the Blender UI property and updates
        the computed AtomGroup if the selection string has changed.

        Notes
        -----
        This is called by update handlers when the user modifies the selection
        through Blender's UI. It ensures the Python-side state reflects
        the PropertyGroup state.
        """
        self.invalidate_ui_cache()
        # Recompute AtomGroup if string changed
        if self.ui_item.string != self._current_string:
            self.apply_selection_string(self.ui_item.string)
            self.sync_to_blender_attribute()

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

    def is_valid(self) -> bool:
        """
        Check if the selection is in a valid state.

        Verifies that all required components are present: the AtomGroup,
        the mask array, and the UI property.

        Returns
        -------
        bool
            True if the selection is valid and usable, False otherwise.

        Examples
        --------
        >>> if not sel.is_valid():
        ...     print("Selection needs reinitialization")
        """
        try:
            has_ag = self._ag is not None
            has_mask = self.mask_array is not None
            has_ui = self.name in self.trajectory.selections.items
            return has_ag and has_mask and has_ui
        except Exception:
            return False

    def ensure_valid(self) -> None:
        """
        Raise an exception if the selection is in an invalid state.

        Useful for catching bugs early in development or providing
        clear error messages when a selection becomes corrupted.

        Raises
        ------
        SelectionError
            If the selection is missing required components or is otherwise invalid.

        Examples
        --------
        >>> sel.ensure_valid()  # Raises if invalid
        >>> # Continue with valid selection operations
        """
        if not self.is_valid():
            has_ag = self._ag is not None
            has_mask = self.mask_array is not None
            try:
                has_ui = self.name in self.trajectory.selections.items
            except Exception:
                has_ui = False

            raise SelectionError(
                f"Selection '{self.name}' is in an invalid state. "
                f"Has AtomGroup: {has_ag}, "
                f"Has mask: {has_mask}, "
                f"In UI: {has_ui}"
            )

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
        - Creates the Python wrapper FIRST and registers it BEFORE creating the UI
          property, to avoid race conditions with update handlers.
        """
        # Extract metadata from the AtomGroup using the class directly
        meta = AtomGroupMetadata(atomgroup)
        name = meta.fallback_name if name == "" else name

        # Create Python wrapper FIRST
        sel = cls(trajectory=trajectory, name=name)
        sel._uuid = meta.uuid
        sel._ag = atomgroup
        sel.mask_array = sel._ag_to_mask()
        sel._current_string = meta.string

        # Register with manager BEFORE creating UI property
        # This prevents race conditions when the UI property setter triggers update handlers
        trajectory.selections.append(sel)

        # Create UI property LAST
        # Note: Setting prop.name may trigger update handlers, but the selection
        # is already registered in _selections, so get() will find it
        prop = trajectory.selections.items.add()
        prop.name = name
        prop.uuid = meta.uuid
        prop.string = meta.string
        prop.updating = meta.updating
        prop.periodic = meta.periodic
        prop.immutable = True

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

    def _validate_name_unique(self, name: str) -> None:
        """
        Validate that a selection name is unique.

        Parameters
        ----------
        name : str
            The name to validate.

        Raises
        ------
        ValueError
            If a selection with this name already exists.
        """
        if name in self._selections:
            raise ValueError(
                f"Selection '{name}' already exists. "
                "Use a unique name or remove the existing selection first."
            )

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

        Raises
        ------
        ValueError
            If a selection with the given name already exists.
        """
        # Note: Not validating uniqueness here for backward compatibility
        # as existing code may rely on overwriting selections
        sel = self._selections[name] = Selection(self.trajectory, name=name)
        sel.add_selection_property(string=string, updating=updating, periodic=periodic)
        return sel

    def create_selection(
        self,
        string: str,
        name: str = "selection_0",
        updating: bool = True,
        periodic: bool = False,
        validate_unique: bool = True,
    ) -> Selection:
        """
        Create and register a new selection with validation.

        This is the preferred method for creating new selections as it
        includes proper validation.

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
        validate_unique : bool, optional
            Whether to validate that the name is unique. Default is True.

        Returns
        -------
        Selection
            The newly created and registered Selection object.

        Raises
        ------
        ValueError
            If validate_unique=True and a selection with the given name already exists.
        """
        if validate_unique:
            self._validate_name_unique(name)

        return self.add(string=string, name=name, updating=updating, periodic=periodic)

    def register_selection(self, selection: Selection) -> None:
        """
        Register an existing Selection object with this manager.

        Adds a Selection object to the internal dictionary without creating
        a new UI property (assumes the Selection was created via a classmethod
        that already handled UI property creation).

        Parameters
        ----------
        selection : Selection
            The Selection object to register.

        Notes
        -----
        This method is typically used internally by factory methods like
        from_atomgroup() after they've created both the UI property and
        the Selection object.
        """
        self._selections[selection.name] = selection

    def append(self, selection: Selection) -> None:
        """
        Register an existing Selection object with this manager.

        .. deprecated::
            Use register_selection() instead for clearer intent.

        Parameters
        ----------
        selection : Selection
            The Selection object to register.
        """
        self.register_selection(selection)

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
        # Note: sel is already registered in _selections by Selection.from_atomgroup()
        # Store the mask as a Blender attribute
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

    def get(self, name: str, lazy_init: bool = True) -> Selection | None:
        """
        Retrieve a Selection by name.

        If the selection exists in the UI properties but not in the Python
        runtime (common when loading from saved .blend files), it will be
        lazily initialized from the UI property.

        Parameters
        ----------
        name : str
            Name of the selection to retrieve.
        lazy_init : bool, optional
            If True, attempts to initialize the selection from UI properties
            if it doesn't exist in Python runtime. Default is True.

        Returns
        -------
        Selection
            The Selection object with the given name.

        Raises
        ------
        KeyError
            If no selection with the given name exists in either Python
            runtime or UI properties.

        Notes
        -----
        Lazy initialization ensures that selections persist correctly when
        loading .blend files, where UI properties are restored but Python
        objects need to be reconstructed.
        """
        # Fast path: already in runtime
        if name in self._selections:
            return self._selections.get(name)

        # Lazy initialization path
        if lazy_init:
            try:
                # Check if UI property exists
                ui_item = self.items[name]
                # Initialize from UI property
                selection = self.from_ui_item(ui_item)
                return selection
            except (KeyError, AttributeError):
                # UI property doesn't exist either
                pass

        # Selection not found anywhere
        raise KeyError(
            f"Selection '{name}' not found. "
            f"Available selections: {list(self._selections.keys())}"
        )

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

    def sync_all_from_ui(self) -> None:
        """
        Synchronize all selections from UI properties.

        Iterates through all UI property items and ensures that each has
        a corresponding Python Selection object. This is useful when loading
        from saved .blend files where the UI properties are restored but
        Python objects need to be reconstructed.

        Notes
        -----
        This method is typically called during trajectory initialization
        or when loading from disk. It ensures that the Python runtime state
        matches the persistent UI property state.
        """
        for ui_item in self.items:
            if ui_item.name not in self._selections:
                try:
                    self.from_ui_item(ui_item)
                except Exception as e:
                    print(
                        f"Warning: Failed to initialize selection '{ui_item.name}': {e}"
                    )

    def has_selection(self, name: str) -> bool:
        """
        Check if a selection exists.

        Checks both the Python runtime and UI properties for the selection.

        Parameters
        ----------
        name : str
            Name of the selection to check.

        Returns
        -------
        bool
            True if the selection exists, False otherwise.
        """
        if name in self._selections:
            return True
        try:
            return name in self.items
        except Exception:
            return False
