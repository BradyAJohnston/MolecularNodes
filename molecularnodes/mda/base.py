"""
Base class for MDAnalysis visualization
"""

from ..session import get_session
from . import utils
from .universes import UniverseManager


class MDAVis:
    """
    Base Class for MDAnalysis visualization

    Attributes
    ----------
    universes: UniverseManager
        Subscriptable and iterable list of UniverseView objects
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
        self.universes = UniverseManager(self.session)
        utils.subscribe_to_active_object()  # subscribe to active object changes
