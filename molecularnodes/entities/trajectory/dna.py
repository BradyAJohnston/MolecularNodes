from MDAnalysis import Universe
import numpy as np
import bpy
from ... import color
from ...blender import coll, nodes
from ... import bpyd
from ...bpyd import AttributeTypes

from .oxdna.OXDNAParser import OXDNAParser
from .oxdna.OXDNAReader import OXDNAReader
from .trajectory import Trajectory

bpy.types.Scene.MN_import_oxdna_topology = bpy.props.StringProperty(
    name="Toplogy",
    description="File path for the topology to import (.top)",
    subtype="FILE_PATH",
    maxlen=0,
)
bpy.types.Scene.MN_import_oxdna_trajectory = bpy.props.StringProperty(
    name="Trajectory",
    description="File path for the trajectory to import (.oxdna / .dat)",
    subtype="FILE_PATH",
    maxlen=0,
)
bpy.types.Scene.MN_import_oxdna_name = bpy.props.StringProperty(
    name="Name",
    description="Name of the created object.",
    default="NewOrigami",
    maxlen=0,
)


class OXDNATrajectory(Trajectory):
    def __init__(self, universe: Universe, world_scale: float = 0.01):
        super().__init__(universe=universe, world_scale=world_scale)
        self._att_names = (
            "base_vector",
            "base_normal",
            "velocity",
            "angular_velocity",
        )

    def create_object(
        self,
        style: str = "vdw",
        name: str = "NewUniverseObject",
        subframes: int = 0,
    ):
        bob = bpyd.create_bob(
            name=name,
            collection=coll.mn(),
            vertices=self.univ_positions * self.world_scale,
            edges=self.bonds,
        )
        self.object = bob.object
        self._update_timestep_values()

        for name in ("chain_id", "res_num", "res_id"):
            self.store_named_attribute(getattr(self, name), name)

    def _update_positions(self, frame: int) -> None:
        super()._update_positions(frame)
        self._update_timestep_values()

    def _update_timestep_values(self):
        for name in self._att_names:
            self.store_named_attribute(
                self.universe.trajectory.ts.data[name] * self.world_scale, name=name
            )


def load(top, traj, name="oxDNA", setup_nodes=True, world_scale=0.01):
    scale_dna = world_scale * 10

    univ = Universe(top, traj, topology_format=OXDNAParser, format=OXDNAReader)
    traj = OXDNATrajectory(univ, world_scale=scale_dna)
    traj.create_object(name=name)
    traj.store_named_attribute(
        data=color.color_chains_equidistant(traj.chain_id),
        name="Color",
        atype=AttributeTypes.FLOAT_COLOR,
    )

    if setup_nodes:
        nodes.create_starting_node_tree(traj.object, style="oxdna", color=None)

    return traj


class MN_OT_Import_OxDNA_Trajectory(bpy.types.Operator):
    bl_idname = "mn.import_oxdna"
    bl_label = "Load"
    bl_description = "Will import the given file and toplogy."
    bl_options = {"REGISTER"}

    def execute(self, context):
        s = context.scene
        load(
            top=s.MN_import_oxdna_topology,
            traj=s.MN_import_oxdna_trajectory,
            name=s.MN_import_oxdna_name,
        )
        return {"FINISHED"}


def panel(layout, scene):
    layout.label(text="Load oxDNA File", icon="FILE_TICK")
    layout.separator()
    row = layout.row()
    row.prop(scene, "MN_import_oxdna_name")
    row.operator("mn.import_oxdna")
    col = layout.column(align=True)
    col.prop(scene, "MN_import_oxdna_topology")
    col.prop(scene, "MN_import_oxdna_trajectory")


# def base_to_int(bases: np.array) -> np.array:
#     """
#     Convert an array of DNA bases to their corresponding MN integer values.

#     Parameters
#     ----------
#     bases : np.array
#         Array of DNA bases.

#     Returns
#     -------
#     np.array
#         Array of corresponding integer values for the DNA bases.
#     """
#     # Values for internal Molecular Nodes use. Defined in data.py
#     base_lookup = {"A": 30, "C": 31, "G": 32, "T": 33}

#     ints = np.array([base_lookup.get(base, -1) for base in bases])

#     return ints


# def is_new_topology(filepath):
#     with open(filepath) as f:
#         firstline = f.readline()

#     return "5 -> 3" in firstline


# def read_topology_new(filepath):
#     with open(filepath, "r") as file:
#         contents = file.read()

#     lines = np.array(contents.split("\n"))

#     def read_seq_line(line):
#         sequence = line.split(" ")[0]
#         return np.array([c for c in sequence])

#     strands = []
#     counter = 0

#     for i, line in enumerate(lines[1:]):
#         bases = read_seq_line(line)
#         arr = np.zeros((len(bases), 4), dtype=int)
#         idx = np.array(range(len(bases)), dtype=int)
#         arr[:, 0] = i + 1  # strand ID
#         arr[:, 1] = base_to_int(bases)  # base
#         bond_3 = idx - 1 + counter
#         bond_5 = idx + 1 + counter
#         bond_3[0] = -1
#         bond_5[-1] = -1
#         arr[:, 2] = bond_3
#         arr[:, 3] = bond_5

#         strands.append(arr)
#         counter += len(bases)

#     return np.vstack(strands)


# def read_topology_old(filepath):
#     """
#     Read the topology from a file and convert it to a numpy array.


#     Strand assignment
#     |  Base assignment
#     |  |  3' Bonded base to the current base (index based on row)
#     |  |  |   5' Bonded base to the current base (index based on row)
#     |  |  |   |
#     S  B  3'  5'
#     S  B  3'  5'
#     S  B  3'  5'

#     Parameters
#     ----------
#     filepath : str
#         The path to the file containing the topology.

#     Returns
#     -------
#     numpy.ndarray
#         The topology as a integer numpy array. Base assignment is (30, 31, 32, 33) where
#         this corresponds to (A, C, G, T) for use inside of Molecular Nodes.

#     """

#     with open(filepath, "r") as file:
#         contents = file.read()

#     lines = np.array(contents.split("\n"))
#     # metadata = lines[0]

#     # read the topology from the file sans the first metadata line
#     # have to initially read as strings, then convert bases to numeric later
#     array_str = np.loadtxt(lines[1:], dtype=str)

#     # convert the columns to numeric
#     array_int = np.zeros(array_str.shape, dtype=int)
#     array_int[:, (0, 2, 3)] = array_str[:, (0, 2, 3)].astype(
#         int
#     )  # easy convert numeric columns to int
#     # convert bases (A, C, G, T) to (30, 31, 32, 33)
#     array_int[:, 1] = base_to_int(array_str[:, 1])

#     return array_int


# def read_trajectory(filepath):
#     """
#     Read an oxDNA trajectory file and return an array of frames.

#     Each frame becomes a 2D array in a stack. Each frame has 5 three-component vectors.
#     The vectors are: (position, base_vector, base_normal, veclocity, angular_velocity),
#     which totals 15 columns in the array. The (velocity, angular_velocity) are optional
#     and can sometimes not appear in the trajectory.

#     Parameters
#     ----------
#     filepath : str
#         The path to the trajectory file.

#     Returns
#     -------
#     frames : ndarray
#         An array of frames, where each frame is a 2D array of positions

#     """
#     # Open the file and read its contents
#     with open(filepath, "r") as file:
#         contents = file.read()

#     # Split the contents into lines
#     lines = np.array(contents.split("\n"))
#     is_meta = np.char.find(lines, "=") > 0

#     group_id = np.cumsum(np.append([True], np.diff(is_meta)))
#     groups = np.unique(group_id)

#     frames = []

#     for group in groups:
#         mask = group == group_id
#         if "=" in lines[mask][0]:
#             continue

#         arr = np.loadtxt(lines[mask])
#         frames.append(arr)

#     return np.stack(frames)


# def store_named_attributes_to_dna_mol(obj, frame, scale_dna=0.1):
#     attributes = ("base_vector", "base_normal", "velocity", "angular_velocity")
#     for i, att in enumerate(attributes):
#         col_idx = np.array([3, 4, 5]) + i * 3

#         try:
#             data = frame[:, col_idx]
#         except IndexError as e:
#             print(f"Unable to get {att} attribute from coordinates. Error: {e}")
#             continue

#         if att != "angular_velocity":
#             data *= scale_dna

#         bpyd.store_named_attribute(
#             obj=obj, data=data, name=att, atype=AttributeTypes.FLOAT_VECTOR
#         )


# def toplogy_to_bond_idx_pairs(topology: np.ndarray):
#     """
#     Convert the given topology array into pairs of indices representing each distinct bond.

#     Strand assignment
#     |  Base assignment
#     |  |  3' Bonded base to the current base (index based on row)
#     |  |  |  5' Bonded base to the current base (index based on row)
#     |  |  |  |
#     1  A -1  1
#     1  G  0  2
#     1  C  1 -1

#     The topology above becomes:
#     np.array([[0, 1], [2, 1]])

#     The order of the bond indices doesn't matter to Blender.

#     Parameters:
#     topology (np.ndarray): Numeric numpy array representing the topology.

#     Returns:
#     np.ndarray: Array of pairs of indices representing each distinct bond.
#     """

#     # to get pairs of indices which represent each distinct bond, which are needed for
#     # edge creation in Blender, take each bonded column and create a 'bond' with itself
#     idx = np.array(list(range(topology.shape[0])))
#     bond_3 = np.vstack((idx, topology[:, 2])).reshape((len(idx), 2))
#     bond_5 = np.vstack((idx, topology[:, 3])).reshape((len(idx), 2))
#     bonds = np.vstack((bond_3, bond_5))

#     # drop where either bond is -1 (not bonded) from the bond indices
#     mask = bonds == -1
#     mask = np.logical_not(mask.any(axis=1))

#     bond_idxs = np.unique(bonds[mask, :], axis=0)

#     return np.sort(bond_idxs, axis=1)


# def load(top, traj, name="oxDNA", setup_nodes=True, world_scale=0.01):
#     # the scale of the oxDNA files seems to be based on nanometres rather than angstrongs
#     # like most structural biology files, so currently adjusting the world_scale to
#     # compensate
#     scale_dna = world_scale * 10

#     # read in the topology and trajectory files
#     is_new_top = is_new_topology(top)
#     if is_new_top:
#         topology = read_topology_new(top)
#     else:
#         topology = read_topology_old(top)

#     trajectory = read_trajectory(traj)
#     n_frames = trajectory.shape[0]

#     # creat toplogy object with positions of the first frame, and the bonds from the
#     # topology object
#     obj = bpyd.create_object(
#         name=name,
#         collection=coll.mn(),
#         vertices=trajectory[0][:, 0:3] * scale_dna,
#         edges=toplogy_to_bond_idx_pairs(topology),
#     )

#     # adding additional toplogy information from the topology and frames objects
#     bpyd.store_named_attribute(
#         obj=obj, data=topology[:, 1], name="res_name", atype=AttributeTypes.INT
#     )
#     bpyd.store_named_attribute(
#         obj=obj, data=topology[:, 0], name="chain_id", atype=AttributeTypes.INT
#     )
#     bpyd.store_named_attribute(
#         obj=obj,
#         data=color.color_chains_equidistant(topology[:, 0]),
#         name="Color",
#         atype=AttributeTypes.FLOAT_COLOR,
#     )
#     store_named_attributes_to_dna_mol(obj, trajectory[0], scale_dna=scale_dna)

#     # if the 'frames' file only contained one timepoint, return the object without creating
#     # any kind of collection for storing multiple frames from a trajectory, and a None
#     # object in place of the frames collection
#     if n_frames == 1:
#         if setup_nodes:
#             nodes.create_starting_node_tree(obj, style="oxdna", color=None)
#         return obj, None

#     # create a collection to store all of the frame objects that are part of the trajectory
#     # they will contain all of the possible attributes which can be interpolated betewen
#     # frames such as position, base_vector, base_normal, velocity, angular_velocity
#     collection = coll.frames(name)
#     for i, frame in enumerate(trajectory):
#         fill_n = int(np.ceil(np.log10(n_frames)))
#         frame_name = f"{name}_frame_{str(i).zfill(fill_n)}"
#         frame_obj = bpyd.create_object(
#             frame[:, 0:3] * scale_dna, name=frame_name, collection=collection
#         )
#         store_named_attributes_to_dna_mol(frame_obj, frame, scale_dna)

#     if setup_nodes:
#         nodes.create_starting_node_tree(
#             obj, coll_frames=collection, style="oxdna", color=None
#         )

#     return obj, collection


# class MN_OT_Import_OxDNA_Trajectory(bpy.types.Operator):
#     bl_idname = "mn.import_oxdna"
#     bl_label = "Load"
#     bl_description = "Will import the given file and toplogy."
#     bl_options = {"REGISTER"}

#     def execute(self, context):
#         s = context.scene
#         load(
#             top=s.MN_import_oxdna_topology,
#             traj=s.MN_import_oxdna_trajectory,
#             name=s.MN_import_oxdna_name,
#         )
#         return {"FINISHED"}


# def panel(layout, scene):
#     layout.label(text="Load oxDNA File", icon="FILE_TICK")
#     layout.separator()
#     row = layout.row()
#     row.prop(scene, "MN_import_oxdna_name")
#     row.operator("mn.import_oxdna")
#     col = layout.column(align=True)
#     col.prop(scene, "MN_import_oxdna_topology")
#     col.prop(scene, "MN_import_oxdna_trajectory")
