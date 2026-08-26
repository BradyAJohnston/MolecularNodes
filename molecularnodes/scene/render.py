import bpy


def enable_gpus(device_type, use_cpus=False):
    preferences = bpy.context.preferences
    cycles_preferences = preferences.addons["cycles"].preferences
    cycles_preferences.refresh_devices()

    activated_gpus = []
    for device in cycles_preferences.devices:
        if device.type == "CPU":
            device.use = use_cpus
        else:
            device.use = device.type == device_type
            if device.use:
                activated_gpus.append(device.name)

    # the compute_device_type enum accepts every backend the build supports,
    # even with no matching hardware present - rendering then aborts, so
    # require an actual device before committing to the backend
    if not activated_gpus:
        raise RuntimeError(f"No {device_type} devices found")

    cycles_preferences.compute_device_type = device_type
    bpy.context.scene.cycles.device = "GPU"

    return activated_gpus


def enable_optimal_gpu():
    options = ["OPTIX", "CUDA", "METAL", "HIP", "ONEAPI"]
    for backend in options:
        try:
            return enable_gpus(backend)
        # TypeError for a backend unsupported by this build, RuntimeError for
        # a supported backend with no devices
        except (TypeError, RuntimeError):
            continue

    raise RuntimeError("Failed to enable GPU backend")
