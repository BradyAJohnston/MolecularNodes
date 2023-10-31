import bpy

def test_menus_registered():
    menu_names = (
        "MN_MT_Node_Color",
        # "MN_MT_Node_Bonds",
        "MN_MT_Node_Style",
        "MN_MT_Node_Select",
        "MN_MT_Node_Assembly",
        "MN_MT_Node_Membranes",
        "MN_MT_Node_DNA",
        "MN_MT_Node_Animate",
        "MN_MT_Node_Utilities", 
        "MN_MT_Node_Density", 
        "MN_MT_Node"
    )
    
    menus = bpy.types.Menu.__subclasses__()
    for name in menu_names:
        assert (f"<class 'molecularnodes.ui.{name}'>" in [str(menu) for menu in menus]) == True
