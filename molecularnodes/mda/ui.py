"""
GUI for MDAnalysis visualization
"""

import bpy


class MDAPanelBase:
    """
    Base mix-in class for all MDAnalysis panels (except Universes)
    """

    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "MDAnalysis"

    @classmethod
    def poll(cls, context):
        """Visible only if a universe is selected in Universes panel"""
        return context.scene.mn.mda.active_index != -1


class MDA_UL_UniversesList(bpy.types.UIList):
    """
    UIList of universes in Universes panel
    """

    def draw_item(
        self,
        context,
        layout,
        data,
        item,
        icon,
        active_data,
        active_property,
        index=0,
        flt_flag=0,
    ):
        custom_icon = "WORLD"
        if self.layout_type in {"DEFAULT", "COMPACT"}:
            row = layout.row()
            name = str(index + 1) + ". "
            split = row.split(factor=0.1)
            col = split.column()
            col.label(text=name)
            col = split.column()
            if item.object is not None:  # can be None if panel is open during draw
                col.prop(item.object, "name", text="", emboss=False)
            hide_icon = "HIDE_OFF" if item.visible else "HIDE_ON"
            row.prop(
                item,
                "visible",
                icon_only=True,
                icon=hide_icon,
            )
            op = row.operator("mda.frame_selected_universe", icon="VIEWZOOM", text="")
            op.index = index
        elif self.layout_type in {"GRID"}:
            layout.alignment = "CENTER"
            layout.label(text="", icon=custom_icon)


class MDA_PT_universes(bpy.types.Panel):
    """
    Panel for list of universes
    """

    bl_idname = "MDA_PT_universe"
    bl_label = "Universes"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "MDAnalysis"

    def draw(self, context):
        layout = self.layout
        sprop = context.scene.mn.mda
        row = layout.row()
        row.template_list(
            "MDA_UL_UniversesList",
            "universes_list",
            sprop,
            "universes",
            sprop,
            "active_index",
            rows=3,
        )
        col = row.column()
        col.operator("mda.add_universe", icon="ADD", text="")
        col.operator("mda.delete_universe", icon="REMOVE", text="")

        row = layout.row()
        row.operator("mda.frame_all", icon="ZOOM_ALL")


class MDA_PT_trajectory(MDAPanelBase, bpy.types.Panel):
    """
    Panel for universe trajectory details
    """

    bl_idname = "MDA_PT_trajectory"
    bl_label = "Trajectory"

    def draw(self, context):
        layout = self.layout
        # To enable the animatate dot next to property in UI
        # layout.use_property_split = True
        # layout.use_property_decorate = True
        object = context.object
        prop = object.mn
        row = layout.row()
        label = "This trajectory has " + str(prop.n_frames) + " frames"
        box = layout.box()
        row.label(text=label)
        row = box.row()
        row.prop(prop, "frame")


CLASSES = [
    MDA_UL_UniversesList,
    MDA_PT_universes,
    MDA_PT_trajectory,
]
