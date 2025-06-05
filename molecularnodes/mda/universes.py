"""
Universes
"""

import bpy
import MDAnalysis as mda
from MDAnalysis.core.groups import AtomGroup
from ..entities.trajectory import Trajectory


class Universe:
    """
    Base Universe class

    Attributes
    ----------
    name: str
        Name of the universe
    object: bpy.types.Object
        Corresponding Blender object

    """

    def __init__(
        self,
        universe: mda.Universe | AtomGroup,
        style: str = "spheres",
        name: str = None,
    ) -> None:
        """Initialize by adding the universe"""
        sprop = bpy.context.scene.mn.mda  # scene pointer property
        uprop = sprop.universes.add()  # add new universe property to collection
        uprop.name = "u" + str(sprop.next_index)  # unique key name
        sprop.next_index += 1
        sprop.active_index = len(sprop.universes) - 1
        object_name = name
        if object_name is None:
            object_name = uprop.name  # use unique name if no name passed
        # create the trajectory
        if isinstance(universe, AtomGroup):
            self.trajectory = Trajectory(universe.universe)
            # TODO: Add style only from the AtomGroup selection
        else:
            self.trajectory = Trajectory(universe)
        # create the blender object
        self.object = self.trajectory.create_object(name=object_name, style=style)
        self.object.mn.update_with_scene = False  # disable by default
        uprop.object = self.object  # link object to property
        oprop = self.object.mn.mda  # object pointer property
        oprop.is_mda_universe = True  # mark as MDA universe
        oprop.universe_key = uprop.name  # save the unique key name for later lookups

    @property
    def name(self) -> str:
        """Get universe name"""
        return self.object.name  # universe name is object name

    @name.setter
    def name(self, value: str) -> None:
        """Set universe name"""
        self.object.name = value

    @property
    def key(self) -> str:
        """Unique key to lookup universes property collection"""
        return self.object.mn.mda.universe_key


class Universes:
    """
    Wrapper class for universes dict to make it subscriptable,
    iterable and support key completions in IPython
    """

    def __init__(self, universes: dict) -> None:
        self._universes = universes

    def __iter__(self):
        """To support iteration"""
        return iter(self._universes.values())

    def __len__(self):
        """Number of universes"""
        return len(self._universes)

    def __getitem__(self, key: str) -> Universe:
        """Return Universe object based on name"""
        return next((u for u in self._universes.values() if u.name == key))

    def _ipython_key_completions_(self):
        """Return dynamic universe (object) names"""
        return [u.name for u in self._universes.values()]
