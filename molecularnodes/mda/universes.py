"""
Classes related to Universe view and management
"""

import bpy
import MDAnalysis as mda
from MDAnalysis.core.groups import AtomGroup
from ..entities.trajectory import Trajectory
from ..session import MNSession
from . import utils


class UniverseView:
    """
    Class responsible for display of a MDAnalysis Universe

    Attributes
    ----------
    name: str
        Name of the universe (same as blender object name)
    trajectory: Trajectory
        MN trajectory entity of the universe
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
        # create the trajectory
        if isinstance(universe, AtomGroup):
            self.trajectory = Trajectory(universe.universe)
            # TODO: Add style only from the AtomGroup selection
        else:
            self.trajectory = Trajectory(universe)
        # create the blender object
        self.object = self.trajectory.create_object(name=name, style=style)
        self.object.mn.update_with_scene = False  # disable by default

    @property
    def name(self) -> str:
        """Get universe name"""
        return self.object.name  # universe name is object name

    @name.setter
    def name(self, value: str) -> None:
        """Set universe name"""
        self.object.name = value

    @property
    def _key(self) -> str:
        """Internal unique key to lookup universes property collection"""
        return self.object.mn.mda.universe_key


class UniverseManager:
    """
    Class responsible for managing universe views (add/delete/clear/list/get)
    """

    _instance = None

    def __new__(cls, _):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self, session: MNSession):
        if self._initialized:
            return
        self._universes = {}  # dict to hold Universe objects created
        self._session = session
        self._initialized = True

    def __iter__(self):
        """To support iteration"""
        return iter(self._universes.values())

    def __len__(self):
        """Number of universes"""
        return len(self._universes)

    def __getitem__(self, key: str) -> UniverseView:
        """Return Universe object based on name"""
        return next((u for u in self._universes.values() if u.name == key))

    def _ipython_key_completions_(self):
        """Return dynamic universe (object) names"""
        return [u.name for u in self._universes.values()]

    def add_universe(
        self,
        universe: mda.Universe | AtomGroup,
        style: str = "spheres",
        name: str = None,
    ) -> UniverseView:
        """Add a Universe or AtomGroup view"""
        sprop = bpy.context.scene.mn.mda  # scene pointer property
        uprop = sprop.universes.add()  # add new universe property to collection
        uprop.name = "u" + str(sprop.next_index)  # unique key name
        sprop.next_index += 1
        sprop.active_index = len(sprop.universes) - 1
        object_name = name
        if object_name is None:
            object_name = uprop.name  # use unique name if no name passed
        # create the universe view object
        u = UniverseView(universe, style, object_name)
        # update universe and object properties
        uprop.object = u.object  # link object to property
        oprop = u.object.mn.mda  # object pointer property
        oprop.is_mda_universe = True  # mark as MDA universe
        oprop.universe_key = uprop.name  # save the unique key name for later lookups
        self._universes[u._key] = u  # add to dict
        print("Added universe", u.name)
        return u

    def delete_universe(self, universe: UniverseView | str) -> None:
        """Delete a universe view"""
        if isinstance(universe, str):
            # lookup universe object based on name
            u = self.get_universe(universe)
            if u is None:
                return
            object = u.object
            key = u._key
        else:
            object = universe.object
            key = universe.object.mn.mda.universe_key
        sprop = bpy.context.scene.mn.mda  # scene pointer property
        index = sprop.universes.find(key)  # find universe property in collection
        if index == -1:
            print(universe, "not found")
            return
        sprop.universes.remove(index)  # remove from universes collection
        sprop.active_index = len(sprop.universes) - 1
        del self._session.entities[object.uuid]  # remove from MNSession
        object_name = object.name  # object / universe name being deleted
        bpy.data.objects.remove(object, do_unlink=True)  # remove Blender object
        del self._universes[key]  # remove universe object from dict
        print("Deleted universe", object_name)

    def clear_universes(self) -> None:
        """Clear all universe views"""
        for u in self._universes.copy().values():
            self.delete_universe(u)

    def list_universes(self) -> None:
        """List all universe views"""
        if len(self._universes) == 0:
            print("No universes")
            return
        for u in self._universes.values():
            print(u.name, u.trajectory)

    def get_universe(self, name: str) -> UniverseView | None:
        """Get universe view based on name"""
        if name not in bpy.data.objects:  # check name in blender objects
            print(name, "universe not found")
            return None
        object = bpy.data.objects[name]
        if not utils.is_mda_universe(object):  # check object is mda universe
            print(name, "is not a mda universe")
            return None
        return self._universes[object.mn.mda.universe_key]
