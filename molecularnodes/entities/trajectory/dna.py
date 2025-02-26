import databpy
from MDAnalysis import Universe

from ... import color
from ...blender import coll, nodes
from ..base import EntityType
from .base import Trajectory

DNA_SCALE = 10


class OXDNA(Trajectory):
    """
    A class to handle oxDNA trajectory data.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        The MDAnalysis Universe object containing the trajectory data
    world_scale : float, optional
        Scaling factor for the world coordinates, by default 0.01

    Attributes
    ----------
    _entity_type : EntityType
        Type of the molecular entity
    _att_names : tuple
        Names of the attributes to track
    """

    def __init__(self, universe: Universe, world_scale: float = 0.01):
        super().__init__(universe=universe, world_scale=world_scale * DNA_SCALE)
        self._entity_type = EntityType.MD_OXDNA
        self._att_names = (
            "base_vector",
            "base_normal",
            "velocity",
            "angular_velocity",
        )

    def _create_object(
        self, style: str | None = "oxdna", name: str = "NewUniverseObject"
    ):
        """
        Create a new object with the trajectory data.

        Parameters
        ----------
        style : str, optional
            Style of the object representation, by default "oxdna"
        name : str, optional
            Name of the new object, by default "NewUniverseObject"

        Returns
        -------
        bpy.types.Object
            The created Blender object
        """
        self.object = databpy.create_object(
            name=name,
            collection=coll.mn(),
            vertices=self.univ_positions,
            edges=self.bonds,
        )
        self._update_timestep_values()

        for name in ("chain_id", "res_id", "res_name"):
            if name == "res_name":
                att_name = "res_num"
            else:
                att_name = name
            self.store_named_attribute(getattr(self, att_name), name)

        self.store_named_attribute(
            data=color.color_chains_equidistant(self.chain_id),
            name="Color",
            atype=databpy.AttributeTypes.FLOAT_COLOR,
        )

        if style:
            nodes.create_starting_node_tree(self.object, style="oxdna", color=None)

        return self.object

    def set_frame(self, frame: int) -> None:
        super()._update_positions(frame)
        self._update_timestep_values()

    def _update_timestep_values(self):
        """
        Update the timestep values for all tracked attributes.
        """
        for name in self._att_names:
            try:
                self.store_named_attribute(
                    self.universe.trajectory.ts.data[name] * self.world_scale, name=name
                )
            except KeyError as e:
                print(e)
