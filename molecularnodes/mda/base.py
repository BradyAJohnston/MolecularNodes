"""
Base class for MDAnalysis visualization
"""

import bpy
import MDAnalysis as mda
from MDAnalysis.core.groups import AtomGroup
from ..session import get_session
from . import utils
from .universes import Universe, Universes


class MDAVis:
    """
    Base Class for MDAnalysis visualization

    Attributes
    ----------
    universes: Universes
        Subscriptable and iterable list of Universe objects
    session: MNSession
        MolecularNodes session object
    """

    def __new__(cls) -> None:
        """This is a singleton class"""
        # try to retrieve object from MNSession
        try:
            session = get_session()
        except AttributeError:
            from .. import register

            register()
            session = get_session()
        if hasattr(session, "MDAVis"):
            return session.MDAVis
        instance = super().__new__(cls)
        session.MDAVis = instance  # add instance to MNSession
        instance._initialized = False
        instance.session = session  # link to MNSession (for deletes)
        return instance

    def __init__(self) -> None:
        if self._initialized:
            return
        self._initialized = True
        self._universes = {}  # dict to hold Universe objects created
        self.universes = Universes(self._universes)
        utils.subscribe_to_active_object()  # subscribe to active object changes

    def add_universe(
        self,
        universe: mda.Universe | AtomGroup,
        style: str = "spheres",
        name: str = None,
    ) -> Universe:
        """Add a Universe or AtomGroup to Blender"""
        u = Universe(universe, style, name)
        self._universes[u.key] = u  # add to dict
        print("Added universe", u.name)
        return u

    def delete_universe(self, universe: Universe | str) -> None:
        """Delete MDA universe from Blender"""
        if isinstance(universe, str):
            # lookup universe object based on name
            u = self.get_universe(universe)
            if u is None:
                return
            object = u.object
            key = u.key
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
        del self.session.entities[object.uuid]  # remove from MNSession
        object_name = object.name  # object / universe name being deleted
        bpy.data.objects.remove(object, do_unlink=True)  # remove Blender object
        del self._universes[key]  # remove universe object from dict
        print("Deleted universe", object_name)

    def list_universes(self) -> None:
        """List all universes"""
        if len(self._universes) == 0:
            print("No universes")
            return
        for u in self._universes.values():
            print(u.name, u.trajectory)

    def get_universe(self, name: str) -> Universe | None:
        """Get universe based on name"""
        if name not in bpy.data.objects:  # check name in blender objects
            print(name, "universe not found")
            return None
        object = bpy.data.objects[name]
        if not utils.is_mda_universe(object):  # check object is mda universe
            print(name, "is not a mda universe")
            return None
        return self._universes[object.mn.mda.universe_key]
