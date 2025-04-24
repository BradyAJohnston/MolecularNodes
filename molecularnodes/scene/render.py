import bpy


def enable_gpus(device_type, use_cpus=False):
    preferences = bpy.context.preferences
    cycles_preferences = preferences.addons["cycles"].preferences
    cycles_preferences.refresh_devices()
    devices = cycles_preferences.devices

    if not devices:
        raise RuntimeError("Unsupported device type")

    activated_gpus = []
    for device in devices:
        if device.type == "CPU":
            device.use = use_cpus
        else:
            device.use = True
            activated_gpus.append(device.name)
            # print("activated gpu", device.name)

    cycles_preferences.compute_device_type = device_type
    bpy.context.scene.cycles.device = "GPU"

    return activated_gpus


def enable_optimal_gpu():
    options = ["OPTIX", "CUDA", "METAL", "HIP", "METAL", "ONEAPI"]
    for backend in options:
        try:
            enable_gpus(backend)
            break
        except TypeError:
            continue

        raise TypeError("No GPU Device found")
