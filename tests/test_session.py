import bpy
import molecularnodes as mn

try:
    mn.unregister()
except Exception:
    pass
mn.register()


def test_session_present():
    assert isinstance(bpy.context.scene.MNSession, mn.session.MNSession)


def test_persistent_handlers_added():
    load_handlers = [handler.__name__ for handler in bpy.app.handlers.load_post]
    save_handlers = [handler.__name__ for handler in bpy.app.handlers.save_post]
    assert "_session_pickle" in save_handlers
    assert "_session_load" in load_handlers
