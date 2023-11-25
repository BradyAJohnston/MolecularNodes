import bpy
import pytest
import molecularnodes as mn

def test_node_name_format():
    assert mn.blender.nodes.format_node_name("MN_style_cartoon") == "Style Cartoon"