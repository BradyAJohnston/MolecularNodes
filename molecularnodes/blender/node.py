from abc import ABCMeta
from typing import Optional, Any
import warnings
import time
import numpy as np
import bpy


class Node(metaclass=ABCMeta):
    def __init__(self, node: bpy.types.Node, chain=[]):

        self.node = node
        self.group = node.id_data
        self.chain = chain

    @property
    def location(self):
        return np.array(self.node.location)

    def new(self, name):
        "Add a new node to the node group."
        try:
            return self.group.nodes.new(f'GeometryNode{name}')
        except RuntimeError:
            return self.group.nodes.new(f'ShaderNode{name}')

    def link(self, name, linkto=0, linkfrom=0):
        "Create a new node along in the chain and create a link to it. Return the new node."
        new_node = self.new(name)
        new_node.location = self.location + np.array((200, 0))

        self.group.links.new(
            self.node.outputs[linkfrom],
            new_node.inputs[linkto]
        )

        return Node(new_node, chain=self.chain + [self])
