"""
Calculations that write attributes on a ``Molecule``

A calculation is a function that takes the ``Universe`` and returns one value
per atom, stored as an attribute under the calculation's name. An updating
calculation runs on every frame change; one that isn't keeps the values from
when it last ran, which is when it was added, refreshed or set to update.

Whether each calculation updates is stored on the object in
``obj.mn.calculations``, so it is saved with the file. The functions live in
Python: ones that can be pickled are restored with a saved session, and
built-in ones such as DSSP are rebuilt from their saved settings.
"""

from __future__ import annotations
import logging
from typing import TYPE_CHECKING, Callable, Iterator
import numpy as np
from databpy import LinkedObjectError
from MDAnalysis import Universe

if TYPE_CHECKING:
    from .base import Molecule

logger = logging.getLogger(__name__)

Calculation = Callable[[Universe], np.ndarray]


class CalculationManager:
    """Calculations that write attributes from the universe.

    Works like a dictionary of functions by attribute name, so
    ``mol.calculations["name"] = func`` adds an updating calculation and
    ``del mol.calculations["name"]`` removes it.

    Parameters
    ----------
    entity : Molecule
        The molecule the calculations write attributes on.
    """

    def __init__(self, entity: Molecule):
        self._entity = entity
        self._funcs: dict[str, Calculation] = {}

    @property
    def ui_items(self):
        """The saved calculation items on the object, ``obj.mn.calculations``."""
        return self._entity.object.mn.calculations

    def _item(self, name: str):
        """The saved item of a calculation, or None without one or an object."""
        try:
            return self.ui_items.get(name)
        except (LinkedObjectError, AttributeError):
            return None

    def add(
        self,
        name: str,
        func: Calculation,
        *,
        updating: bool = True,
        refresh: bool = True,
    ) -> None:
        """Add a calculation, replacing any with the same name.

        A calculation that isn't updating is run once now, an updating one on
        the next frame change.

        Parameters
        ----------
        name : str
            Name of the attribute the calculation writes.
        func : Callable[[Universe], np.ndarray]
            Function returning one value per atom.
        updating : bool, default=True
            Whether to recalculate on every frame change.
        refresh : bool, default=True
            Whether to run a calculation that isn't updating now. Without this
            the attribute keeps its current values until the next refresh.
        """
        self._funcs[name] = func
        try:
            items = self.ui_items
        except (LinkedObjectError, AttributeError):
            # no object to store attributes on yet
            return
        item = items.get(name)
        if item is None:
            item = items.add()
            item.name = name
        if item.updating != updating:
            item.updating = updating
        if refresh and not updating:
            self.refresh(name)
        self._entity._sync_frame_attributes()

    def updating(self, name: str) -> bool:
        """Whether a calculation recalculates on every frame change."""
        item = self._item(name)
        if item is None:
            if name in self._funcs:
                return True
            raise KeyError(name)
        return item.updating

    def set_updating(self, name: str, updating: bool) -> None:
        """Set whether a calculation recalculates on every frame change."""
        item = self._item(name)
        if item is None:
            raise KeyError(name)
        if item.updating != updating:
            item.updating = updating

    def refresh(self, name: str) -> None:
        """Run a calculation now and store its values, whether or not it updates."""
        values = self._values([name])
        self._entity._store_frame_attributes(values, "CALCULATION")

    def _values(self, names) -> dict[str, np.ndarray]:
        """Run the named calculations, logging any that fail."""
        values = {}
        for name in names:
            try:
                values[name] = self._funcs[name](self._entity.universe)
            except Exception as e:
                logger.error(
                    f"Failed to update calculation '{name}': {e}", exc_info=True
                )
        return values

    def _frame_values(self) -> dict[str, np.ndarray]:
        """Values of the updating calculations for the current frame."""
        names = [
            name
            for name in self._funcs
            if (item := self._item(name)) is None or item.updating
        ]
        return self._values(names)

    def __setitem__(self, name: str, func: Calculation) -> None:
        self.add(name, func)

    def __getitem__(self, name: str) -> Calculation:
        return self._funcs[name]

    def __delitem__(self, name: str) -> None:
        del self._funcs[name]
        try:
            items = self.ui_items
        except (LinkedObjectError, AttributeError):
            return
        index = items.find(name)
        if index != -1:
            items.remove(index)
        self._entity._sync_frame_attributes()

    def __contains__(self, name: object) -> bool:
        return name in self._funcs

    def __iter__(self) -> Iterator[str]:
        return iter(self._funcs)

    def __len__(self) -> int:
        return len(self._funcs)

    def items(self):
        return self._funcs.items()

    def get(self, name: str, default=None):
        return self._funcs.get(name, default)
