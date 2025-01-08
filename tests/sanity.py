import bpy


def current_frame() -> int:
    return bpy.context.scene.frame_current


def changing_frame(scene):
    print(f"{current_frame()=}")


bpy.app.handlers.frame_change_pre.append(changing_frame)

print(f"{current_frame()=}")

bpy.context.scene.frame_set(3)

print(f"{current_frame()=}")
