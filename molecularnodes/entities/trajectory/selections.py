import MDAnalysis as mda
import numpy.typing as npt
import numpy as np


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

    @classmethod
    def from_atomgroup(cls, atomgroup: mda.AtomGroup, name: str = ""):
        "Create a Selection object from an AtomGroup"

        # set default value
        selection_str = f"sel_{atomgroup.n_atoms}_atoms"
        updating = False
        periodic = False

        # if class is an UpdatingAtomGroup
        if atomgroup.__class__.__name__ == "UpdatingAtomGroup":
            updating = True
            # assuming it's a single selection
            # MDA do support `u.select_atoms('index 0', 'around 5 index 0')
            selection_str = atomgroup._selection_strings[0]
            periodic = False
            try:
                if atomgroup._selections[0].periodic:
                    periodic = True
            except AttributeError as e:
                print(e)

        if name == "":
            name = selection_str
        selection = cls(atomgroup.universe, "all", name, updating, periodic)

        selection.selection_str = selection_str
        selection.ag = atomgroup
        selection.mask_array = selection._ag_to_mask()
        return selection
