"""
Calculations that write attributes on a ``Molecule``

A calculation is a function that takes the ``Universe`` and returns one value
per atom, stored as an attribute under the calculation's name, or one that
writes a data object itself and returns None, such as hydrogen bonds. An updating
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
import bpy
import numpy as np
from databpy import LinkedObjectError
from MDAnalysis import Universe
from .dssp import DSSPCalculation
from .hbonds import HBondCalculation

if TYPE_CHECKING:
    from .base import Molecule

logger = logging.getLogger(__name__)

Calculation = Callable[[Universe], np.ndarray | None]


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
        self._dssp = DSSPCalculation(entity)
        self._hbonds = HBondCalculation(entity)
        self._failed_texts: set[str] = set()

    def dssp(
        self, selection: str | None = None, *, updating: bool | None = None
    ) -> DSSPCalculation:
        """Secondary structure from DSSP, as the `sec_struct` calculation.

        Sets DSSP up if it isn't already, and changes only the settings given.
        A trajectory whose files have no secondary structure has it set up on
        load, not updating, so it has run on the first frame.

        Parameters
        ----------
        selection : str, optional
            MDAnalysis selection string limiting the atoms DSSP is run on,
            `"protein"` when setting up. Use this to exclude termini or
            non-standard residues that lack the required (N, CA, C, O) atoms.
            Protein residues outside it are drawn as loops.
        updating : bool, optional
            Whether to recalculate on every frame change, `True` when setting
            up. Otherwise the values are kept from when it last ran.

        Returns
        -------
        DSSPCalculation
            The DSSP calculation, whose methods choose what is shown, such as
            `show_sliding_window_average()`

        Examples
        --------
        ```python
        dssp = mol.calculations.dssp(updating=True)
        dssp.show_trajectory_average(threshold=0.5)
        ```
        """
        self._dssp._init(selection=selection, updating=updating)
        return self._dssp

    def hbonds(
        self, selection: str | None = None, *, updating: bool | None = None, **kwargs
    ) -> HBondCalculation:
        """Hydrogen bonds from MDAnalysis, as the `hbonds` calculation.

        Sets hydrogen bonds up if they aren't already, and changes only the
        settings given. The calculation writes the bonds as edges on a data
        object, `HBondCalculation.object`, for the `Style HBonds` node to draw.

        Parameters
        ----------
        selection : str, optional
            MDAnalysis selection string limiting the atoms considered as
            hydrogens, donors and acceptors, `"all"` when setting up.
        updating : bool, optional
            Whether to recalculate on every frame change, `True` when setting
            up. Otherwise the bonds are kept from when it last ran.
        **kwargs
            Passed to `MDAnalysis.analysis.hydrogenbonds.HydrogenBondAnalysis`,
            overriding the defaults and any given before. The cutoffs
            `d_a_cutoff` (3.0 Å), `d_h_a_angle_cutoff` (150 degrees) and
            `d_h_cutoff` (1.2 Å) are stored as properties, adjustable in the
            Attributes panel, and see `HBondCalculation.cutoffs`. By default hydrogens
            (H or D), donors and acceptors (N or O) are selected by element,
            donors are paired with hydrogens by distance, and selections are not
            updated per frame.

        Returns
        -------
        HBondCalculation
            The hydrogen bond calculation, whose `object` is the `Data Object`
            for a `Style HBonds` node

        Raises
        ------
        ValueError
            When there are no hydrogens, which the analysis requires

        Examples
        --------
        ```python
        hbonds = mol.calculations.hbonds(d_a_cutoff=3.5)
        with mol.tree.reset() as (atoms, join):
            atoms >> StyleHBonds(data_object=hbonds.object) >> join
        ```
        """
        self._hbonds._init(selection=selection, updating=updating, **kwargs)
        return self._hbonds

    def from_text(
        self, text: bpy.types.Text | str, *, updating: bool = True
    ) -> list[str]:
        """Add the `BaseCalculation` subclasses defined in a text block.

        The text is run as Python, with `BaseCalculation` and `np` already
        available, and each subclass it defines is added as a calculation
        named after the subclass, replacing any with the same name.

        Parameters
        ----------
        text : bpy.types.Text or str
            The text block, or its name in `bpy.data.texts`.
        updating : bool, default=True
            Whether the calculations recalculate on every frame change.

        Returns
        -------
        list[str]
            The names of the calculations that were added.

        Raises
        ------
        ValueError
            When the text defines no `BaseCalculation` subclasses.
        """
        if isinstance(text, str):
            text = bpy.data.texts[text]
        calculations = calculations_from_code(text.as_string(), self._entity)
        if not calculations:
            raise ValueError(f"No BaseCalculation subclasses found in '{text.name}'")
        n_atoms = self._entity.atoms.n_atoms
        values = {}
        for calculation in calculations:
            result = calculation.calculate()
            if len(result) != n_atoms:
                raise ValueError(
                    f"'{calculation.name}' returned {len(result)} values, "
                    f"expected one for each of the {n_atoms} atoms"
                )
            values[calculation.name] = result
        for calculation in calculations:
            self.add(calculation.name, calculation, updating=updating, refresh=False)
            item = self._item(calculation.name)
            if item is not None:
                item.text = text
        self._entity._store_frame_attributes(values, "CALCULATION")
        return list(values)

    def _restore_texts(self) -> None:
        """Load calculations from their text blocks again after a file is reopened.

        Calculations loaded from text can't be saved with the file, but their
        items on the object record the text block, so they are loaded from it
        again with the updating flags they had.
        """
        try:
            items = self.ui_items
        except (LinkedObjectError, AttributeError):
            return
        missing = [
            item
            for item in items
            if item.text is not None and item.name not in self._funcs
        ]
        for text in {item.text for item in missing}:
            if text.name in self._failed_texts:
                continue
            updating = {item.name: item.updating for item in missing}
            try:
                names = self.from_text(text)
            except Exception as e:
                # don't retry on every frame change
                self._failed_texts.add(text.name)
                logger.warning(f"Could not load calculations from '{text.name}': {e}")
                continue
            for name in names:
                self.set_updating(name, updating.get(name, True))

    @property
    def ui_items(self):
        """The saved calculation items on the object, ``obj.mn.calculations``."""
        return self._entity.props.calculations

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
        data_object: bpy.types.Object | None = None,
    ) -> None:
        """Add a calculation, replacing any with the same name.

        A calculation that isn't updating is run once now, an updating one on
        the next frame change.

        Parameters
        ----------
        name : str
            Name of the attribute the calculation writes.
        func : Callable[[Universe], np.ndarray | None]
            Function returning one value per atom, or writing `data_object`
            and returning None.
        updating : bool, default=True
            Whether to recalculate on every frame change.
        refresh : bool, default=True
            Whether to run a calculation that isn't updating now. Without this
            the attribute keeps its current values until the next refresh.
        data_object : bpy.types.Object, optional
            Object the calculation writes rather than an attribute, recorded
            with its settings on the object.
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
        if item.object != data_object:
            item.object = data_object
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
                value = self._funcs[name](self._entity.universe)
            except Exception as e:
                logger.error(
                    f"Failed to update calculation '{name}': {e}", exc_info=True
                )
                continue
            # calculations that write an object return no attribute values
            if value is not None:
                values[name] = value
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
        if name == "sec_struct":
            # removing DSSP's calculation removes DSSP, it can be set up again
            self._dssp = DSSPCalculation(self._entity)
        elif name == "hbonds":
            # its data object keeps the bonds it last had
            self._hbonds = HBondCalculation(self._entity)
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


class BaseCalculation:
    """
    A base class for custom calculations that are run on every frame.

    Create a subclass and write the calculation inside its `calculate()`
    method. The result is stored as an attribute on the molecule, named
    after the subclass, with its type inferred from the returned array.
    """

    def __init__(self, molecule: Molecule, name: str | None = None):
        self.molecule = molecule
        self.name = name or type(self).__name__

    @property
    def universe(self) -> Universe:
        return self.molecule.universe

    def calculate(self) -> np.ndarray:
        """
        Calculate the values for the current frame.

        Use `self.universe` to access the data for the current frame.

        Returns
        -------
        np.ndarray
            One value for each atom, so the length must match the number
            of atoms in the molecule.
        """
        raise NotImplementedError("Subclasses must implement calculate()")

    def __call__(self, universe: Universe) -> np.ndarray:
        return self.calculate()


def calculations_from_code(code: str, molecule: Molecule) -> list[BaseCalculation]:
    namespace = {"BaseCalculation": BaseCalculation, "np": np}
    exec(code, namespace)
    return [
        obj(molecule)
        for obj in namespace.values()
        if isinstance(obj, type)
        and issubclass(obj, BaseCalculation)
        and obj is not BaseCalculation
    ]
