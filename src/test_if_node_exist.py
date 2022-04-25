import bpy

new_node = new_node
need_to_append = False
try: 
    bpy.data.node_groups.get(new_node).name == new_node
except:
    need_to_append = True