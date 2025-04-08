from uuid import uuid1
import MDAnalysis as mda
import numpy as np
import numpy.typing as npt


class Selection:
    def __init__(self, trajectory, name: str = "selection_0"):
        self._ag: mda.AtomGroup | None = None
        self._uuid: str = str(uuid1())
        self._name = name
        self._current_selection_str: str = ""
        self.trajectory = trajectory

    def add_selection_property(
        self, selection_str: str = "all", updating=True, periodic=True
    ):
        prop = self.trajectory.object.mn_trajectory_selections.add()
        prop.name = self.name
        prop.uuid = self._uuid
        self.updating = updating
        self.periodic = periodic
        self.set_atom_group(selection_str)
        self.selection_str = selection_str
        self.mask_array = self._ag_to_mask()

    @property
    def name(self) -> str:
        return self._name

    @property
    def ui_item(self):
        return self.trajectory.object.mn_trajectory_selections[self.name]

    @property
    def selection_str(self) -> str:
        return self.ui_item.selection_str

    @selection_str.setter
    def selection_str(self, selection_str: str) -> None:
        self.ui_item.selection_str = selection_str
        try:
            self.set_atom_group(selection_str)
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

    def set_atom_group(self, selection_str: str) -> None:
        if self._current_selection_str == selection_str:
            return
        try:
            self._ag = self.trajectory.universe.select_atoms(
                selection_str, updating=self.updating, periodic=self.periodic
            )
            self.mask_array = self._ag_to_mask()
            self._current_selection_str = selection_str
            self.message = ""
        except Exception as e:
            self._current_selection_str = selection_str
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
        self.trajectory.set_boolean(self.to_mask(), name=self.name)

    def to_mask(self) -> npt.NDArray[np.bool_]:
        "Returns the selection as a 1D numpy boolean mask. If updating=True, recomputes selection."
        if self.updating:
            self.mask_array = self._ag_to_mask()
        return self.mask_array

    @classmethod
    def from_atomgroup(cls, trajectory, atomgroup: mda.AtomGroup, name: str = ""):
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
            try:
                if atomgroup._selections[0].periodic:
                    periodic = True
            except AttributeError:
                # some selections don't have the periodic attribute
                pass
            except Exception as e:
                print(e)

        if name == "":
            name = selection_str
        selection = cls(trajectory=trajectory, name=name)
        trajectory.selections[selection.name] = selection

        prop = trajectory.object.mn_trajectory_selections.add()
        prop.name = name
        prop.uuid = selection._uuid

        selection._ag = atomgroup
        selection.mask_array = selection._ag_to_mask()
        selection._current_selection_str = name
        selection.updating = updating
        selection.periodic = periodic
        selection.immutable = True
        selection.selection_str = name

        return selection
