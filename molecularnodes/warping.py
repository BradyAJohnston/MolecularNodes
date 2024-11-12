import bpy
import warp as wp
import numpy as np
from . import bpyd


class WarpSimulator:
    def __init__(self, num_particles):
        # Initialize Warp
        wp.init()

        # Initialize random positions and velocities on CPU first
        positions_np = np.random.randn(num_particles, 3).astype(np.float32)
        velocities_np = np.zeros((num_particles, 3), dtype=np.float32)

        # Create simulation parameters
        self.num_particles = num_particles

        # Create Warp arrays from NumPy arrays
        self.positions = wp.from_numpy(positions_np, dtype=wp.vec3, device="cuda")
        self.velocities = wp.from_numpy(velocities_np, dtype=wp.vec3, device="cuda")

        # Create mesh object for visualization
        self.create_particle_mesh()

        # Initialize simulation parameters
        self.dt = 0.01  # timestep

    @wp.kernel
    def simulate_step(
        positions: wp.array(dtype=wp.vec3),
        velocities: wp.array(dtype=wp.vec3),
        dt: float,
    ):
        tid = wp.tid()

        g = wp.vec3(0.0, -9.81, 0.0)
        velocities[tid] = velocities[tid] + g * dt
        positions[tid] = positions[tid] + velocities[tid] * dt

    def create_particle_mesh(self):
        # Create mesh for particles
        name = "ParticleObject"
        try:
            self.bob = bpyd.BlenderObject(bpy.data.objects[name])
        except KeyError:
            self.bob = bpyd.create_bob(
                np.zeros((self.num_particles, 3), float), name=name
            )

    def step_simulation(self):
        # Run simulation step on GPU
        wp.launch(
            kernel=self.simulate_step,
            dim=self.num_particles,
            inputs=[self.positions, self.velocities, self.dt],
        )

        # Copy positions back to CPU
        positions_np = self.positions.numpy()

        self.bob.position = positions_np
        self.bob.store_named_attribute(self.velocities.numpy(), "velocity")


def frame_change_handler(scene):
    if hasattr(bpy.context.scene, "warp_simulator"):
        bpy.context.scene.warp_simulator.step_simulation()


def register():
    # Create simulator instance
    num_particles = 1000  # Adjust as needed
    simulator = WarpSimulator(num_particles)

    # Store simulator instance in scene
    bpy.types.Scene.warp_simulator = simulator

    # Register frame change handler
    bpy.app.handlers.frame_change_post.append(frame_change_handler)


def unregister():
    if hasattr(bpy.context.Scene, "warp_simulator"):
        del bpy.types.Scene.warp_simulator
    bpy.app.handlers.frame_change_post.remove(frame_change_handler)


if __name__ == "__main__":
    register()
