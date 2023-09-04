import MolecularNodes as mn
import bpy
import pytest

def test_menus_registered():
    
    menu_names = (
        "MN_MT_Add_Node_Menu_Properties",
        "MN_MT_Add_Node_Menu_Color",
        "MN_MT_Add_Node_Menu_Bonds",
        "MN_MT_Add_Node_Menu_Styling",
        "MN_MT_Add_Node_Menu_Selections",
        "MN_MT_Add_Node_Menu_Assembly",
        "MN_MT_Add_Node_Menu_Membranes",
        "MN_MT_Add_Node_Menu_DNA",
        "MN_MT_Add_Node_Menu_Animation",
        "MN_MT_Add_Node_Menu_Utilities", 
        "MN_MT_Add_Node_Menu_Density", 
        "MN_MT_Add_Node_Menu"
    )
    
    menus = bpy.types.Menu.__subclasses__()
    menu_is_present = True
    
    for name in menu_names:
        assert (f"<class 'MolecularNodes.ui.{name}'>" in [str(menu) for menu in menus]) == True
