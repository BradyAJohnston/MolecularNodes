"""Run inside GUI Blender: export every Grease Pencil object in the camera view to SVG."""

import sys
import bpy

out = sys.argv[sys.argv.index("--") + 1]


def run():
    win = bpy.context.window_manager.windows[0]
    area = max(win.screen.areas, key=lambda a: a.width * a.height)
    area.type = "VIEW_3D"
    region = next(r for r in area.regions if r.type == "WINDOW")
    space = area.spaces.active
    space.region_3d.view_perspective = "CAMERA"
    with bpy.context.temp_override(
        window=win, area=area, region=region, space_data=space
    ):
        r = bpy.ops.wm.grease_pencil_export_svg(
            filepath=out,
            selected_object_type="VISIBLE",
            use_fill=True,
            use_clip_camera=True,
            use_uniform_width="--uniform" in sys.argv,
        )
    print("EXPORT", r, out)
    bpy.ops.wm.quit_blender()


bpy.app.timers.register(run, first_interval=1.0)
