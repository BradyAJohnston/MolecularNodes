import MolecularNodes as mn
import bpy
import pytest

def test_menus_registered():
    
    menu_names = (
        "MOL_MT_Add_Node_Menu_Properties",
        "MOL_MT_Add_Node_Menu_Color",
        "MOL_MT_Add_Node_Menu_Bonds",
        "MOL_MT_Add_Node_Menu_Styling",
        "MOL_MT_Add_Node_Menu_Selections",
        "MOL_MT_Add_Node_Menu_Assembly",
        "MOL_MT_Add_Node_Menu_Membranes",
        "MOL_MT_Add_Node_Menu_DNA",
        "MOL_MT_Add_Node_Menu_Animation",
        "MOL_MT_Add_Node_Menu_Utilities", 
        "MOL_MT_Add_Node_Menu_Density", 
        "MOL_MT_Add_Node_Menu"
    )
    
    menus = bpy.types.Menu.__subclasses__()
    menu_is_present = True
    
    for name in menu_names:
        assert (f"<class 'MolecularNodes.ui.{name}'>" in [str(menu) for menu in menus]) == True
