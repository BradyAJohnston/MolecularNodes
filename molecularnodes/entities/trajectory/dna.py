from typing import Set

import bpy
import databpy
from MDAnalysis import Universe

from ... import color
from ...blender import coll, nodes
from ..base import EntityType
from .base import Trajectory
from .ops import TrajectoryImportOperator
from .oxdna.OXDNAParser import OXDNAParser
from .oxdna.OXDNAReader import OXDNAReader

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

    def _create_object(self, style: str = "oxdna", name: str = "NewUniverseObject"):
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


def load(top, traj, name="oxDNA", style="oxdna", world_scale=0.01):
    """
    Load an oxDNA trajectory.

    Parameters
    ----------
    top : str
        Path to topology file
    traj : str
        Path to trajectory file
    name : str, optional
        Name for the created object, by default "oxDNA"
    style : str, optional
        Style of representation, by default "oxdna"
    world_scale : float, optional
        Scaling factor for world coordinates, by default 0.01

    Returns
    -------
    OXDNA
        The created trajectory object
    """
    univ = Universe(top, traj, topology_format=OXDNAParser, format=OXDNAReader)
    traj = OXDNA(univ, world_scale=world_scale)
    traj.create_object(name=name, style=style)
    return traj


class MN_OT_Import_OxDNA_Trajectory(TrajectoryImportOperator):
    """
    Blender operator for importing oxDNA trajectories.
    """

    bl_idname = "mn.import_oxdna"

    def execute(self, context: bpy.types.Context | None) -> Set[str]:
        load(top=self.topology, traj=self.trajectory, name=self.name)
        return {"FINISHED"}


def panel(layout: bpy.types.UILayout, scene: bpy.types.Scene) -> None:
    """
    Create the panel layout for oxDNA import.

    Parameters
    ----------
    layout : bpy.types.UILayout
        Layout to add elements to
    scene : bpy.types.Scene
        Current scene
    """
    layout.label(text="Load oxDNA File", icon="FILE_TICK")
    layout.separator()
    row = layout.row()
    row.prop(scene.mn, "import_oxdna_name")
    op = row.operator("mn.import_oxdna")
    op.name = scene.mn.import_oxdna_name
    op.topology = scene.mn.import_oxdna_topology
    op.trajectory = scene.mn.import_oxdna_trajectory
    col = layout.column(align=True)
    col.prop(scene.mn, "import_oxdna_topology")
    col.prop(scene.mn, "import_oxdna_trajectory")
