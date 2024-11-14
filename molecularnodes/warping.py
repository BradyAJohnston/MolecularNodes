import bpy
import warp as wp
from warp import sim
import numpy as np
from . import bpyd
# from math import sin, cos

from .entities import fetch


class WarpSimulator:
    def __init__(self, num_particles):
        # Initialize Warp
        wp.init()

        # Initialize simulation parameters
        self.num_particles = num_particles
        self.radius = 0.15
        self.frame_dt = 1.0 / 24  # 60 fps
        self.scale = 1

        # Create builder for simulation
        builder = sim.ModelBuilder(up_vector=wp.vec3(0, 0, 1))
        builder.default_particle_radius = self.radius

        n_x = int(num_particles ** (1 / 3))
        n_y = int(num_particles ** (1 / 3))
        n_z = int(num_particles ** (1 / 3))

        b = builder.add_body()

        builder.add_shape_box(
            pos=wp.vec3(0, 0, 0),
            hx=1 * self.scale,
            hy=1 * self.scale,
            hz=1 * self.scale,
            body=b,
        )

        self.mol = fetch("4OZS", ca_only=True, style="ribbon")

        for pos in self.mol.position:
            builder.add_particle(
                pos=pos, vel=wp.vec3(0, 0, 0), radius=self.radius, mass=1
            )

        self.bob = self.mol.bob

        # for edge in self.mol.edges:
        #     i, j = edge.vertices
        #     builder.add_spring(i, j, 1e2, 0.0, 0.0)
        edges = []
        chain_id = self.mol.array.chain_id
        array = self.mol.position
        for i, pos1 in enumerate(array):
            for j, pos2 in enumerate(array):
                if chain_id[i] != chain_id[j]:
                    continue
                if j <= i:
                    continue

                distance = np.linalg.norm(pos1 - pos2)
                if distance > 2:
                    continue
                edges.append([i, j])
                builder.add_spring(i, j, 10, 0.0, 0.0)

        # visualise the springs as edges in the mesh
        positions = self.bob.position.copy()
        # self.bob.object.data.clear_geometry()
        # self.bob.object.data.from_pydata(positions, edges, [])
        # for

        # builder.add_particle_grid(
        #     dim_x=n_x,
        #     dim_y=n_y,
        #     dim_z=n_z,
        #     cell_x=0.1 * 2.0,
        #     cell_y=0.1 * 2.0,
        #     cell_z=0.1 * 2.0,
        #     pos=wp.vec3(-1.0, 0.0, 0.0),
        #     rot=wp.quat_identity(),
        #     vel=wp.vec3(0.0, 0.0, 10.0),
        #     mass=1,
        #     jitter=self.radius * 0.1,
        # )

        # self.num_particles = n_x * n_y * n_z

        # for i in range(self.num_particles):
        #     builder.add_spring(i - 1, i, 3, 0.0, 0)

        # # Create particles in a spiral
        # for i in range(num_particles):
        #     builder.add_particle(
        #         pos=wp.vec3(sin(i / 10) * 3, cos(i / 10) * 3, i / 10),
        #         vel=wp.vec3(0, 0, 0),
        #         radius=0.1,
        #         mass=1,
        #     )
        # builder.add_spring(i - 1, i, 1.0e2, 0.0, 0)

        # Finalize and build the model
        self.model = builder.finalize("cuda")
        # Create states
        self.state_0 = self.model.state()
        self.state_1 = self.model.state()

        # Create integrator
        self.integrator = wp.sim.XPBDIntegrator(10)

        # Create mesh object for visualization
        # self.create_particle_mesh()

    def create_particle_mesh(self):
        name = "ParticleObject"
        try:
            self.bob = bpyd.BlenderObject(bpy.data.objects[name])
            self.bob.position = self.position
        except KeyError:
            self.bob = bpyd.create_bob(self.position, name=name)

    @property
    def position(self) -> np.ndarray:
        return self.state_0.particle_q.numpy()

    @property
    def velocity(self) -> np.ndarray:
        return self.state_0.particle_qd.numpy()

    def simulate(self):
        self.state_0.clear_forces()
        self.state_1.clear_forces()
        self.model.particle_grid.build(self.state_0.particle_q, self.radius * 2)
        wp.sim.collide(self.model, self.state_0)
        self.integrator.simulate(self.model, self.state_0, self.state_1, 1 / 24)
        # swap states
        (self.state_0, self.state_1) = (self.state_1, self.state_0)

    def object_as_wp_transform(self, obj: bpy.types.Object):
        return wp.transform(wp.vec3(obj.location), wp.quat(obj.rotation_quaternion))

    def step(self):
        # get the blender object the user is interacting with and extract
        # the transformations from it
        obj = bpy.data.objects["Cube"]
        if obj.rotation_mode == "QUATERNION":
            rot = obj.rotation_quaternion
        elif obj.rotation_mode == "XYZ":
            rot = obj.rotation_euler.to_quaternion()
        else:
            raise ValueError(f"Unsupported rotation {obj.rotation_type}")

        transform = wp.transform(wp.vec3(*obj.location), wp.quat(*rot))
        self.state_0.body_q.assign([transform])
        # Run simulation substeps
        self.simulate()

        # Update Blender mesh
        self.bob.position = self.position
        self.bob.store_named_attribute(self.velocity, "velocity")


def frame_change_handler(scene):
    if hasattr(bpy.context.scene, "warp_simulator"):
        bpy.context.scene.warp_simulator.step()


def register():
    # Create simulator instance
    num_particles = 10_000  # Adjust as needed
    simulator = WarpSimulator(num_particles)

    # Store simulator instance in scene
    bpy.types.Scene.warp_simulator = simulator

    # Register frame change handler
    bpy.app.handlers.frame_change_post.append(frame_change_handler)


def unregister():
    if hasattr(bpy.context.Scene, "warp_simulator"):
        del bpy.types.Scene.warp_simulator
    bpy.app.handlers.frame_change_post.remove(frame_change_handler)
