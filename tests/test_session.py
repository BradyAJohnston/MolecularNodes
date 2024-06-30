import bpy
import molecularnodes as mn

mn._test_register()


def test_session_present():
    assert isinstance(bpy.context.scene.MNSession, mn.session.MNSession)


def test_persistent_handlers_added():
    load_handlers = [handler.__name__ for handler in bpy.app.handlers.load_post]
    save_handlers = [handler.__name__ for handler in bpy.app.handlers.save_post]
    assert "_pickle" in save_handlers
    assert "_load" in load_handlers
