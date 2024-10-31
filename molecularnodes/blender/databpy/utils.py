import bpy


def evaluate_object(obj: bpy.types.Object):
    "Return an object which has the modifiers evaluated."
    obj.update_tag()
    return obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
