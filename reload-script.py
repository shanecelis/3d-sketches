#!/usr/bin/env python3
# Stick this in Blender's scripting window so you can easily reload/reimport
# from a file. Compatible with Blender 4.5+.

import bpy
from pathlib import Path

# Set this to your .glb path if auto-detection fails (e.g. when running from Text editor).
# Leave as None to auto-detect from script location or blend file location.
# GLB_PATH_OVERRIDE = None  # e.g. "/Users/shane/Projects/3d-sketches/assets/assemble-rocket-ghost.glb"
GLB_PATH_OVERRIDE = "/Users/shane/Projects/3d-sketches/assets/assemble-rocket-ghost.glb"

COLLECTION_NAME = "GLB_RELOAD"

def _resolve_glb_path():
    if GLB_PATH_OVERRIDE:
        p = Path(GLB_PATH_OVERRIDE).expanduser().resolve()
        if p.exists():
            return str(p)
        return str(p)  # still pass through so error message shows the path
    try:
        script_dir = Path(__file__).resolve().parent
    except NameError:
        if bpy.data.filepath:
            script_dir = Path(bpy.path.abspath("//"))
        else:
            script_dir = Path.home() / "Projects" / "3d-sketches"
    p = (script_dir / "assets" / "assemble-rocket-ghost.glb").resolve()
    return str(p)

glb_path = _resolve_glb_path()
if not Path(glb_path).exists():
    raise FileNotFoundError(
        f"GLB not found: {glb_path}\n"
        "Save your .blend in the project folder, or set GLB_PATH_OVERRIDE at the top of this script."
    )

# Delete previous import (iterate over list copy to avoid modifying during iteration)
if COLLECTION_NAME in bpy.data.collections:
    col_prev = bpy.data.collections[COLLECTION_NAME]
    for obj in list(col_prev.objects):
        bpy.data.objects.remove(obj, do_unlink=True)

# Import (Blender 4.5 requires a valid absolute path)
bpy.ops.import_scene.gltf(filepath=glb_path)

# Move imported objects into a known collection
col = bpy.data.collections.get(COLLECTION_NAME)
if not col:
    col = bpy.data.collections.new(COLLECTION_NAME)
    bpy.context.scene.collection.children.link(col)

for obj in list(bpy.context.selected_objects):
    if obj.name not in col.objects:
        col.objects.link(obj)
    try:
        bpy.context.scene.collection.objects.unlink(obj)
    except RuntimeError:
        pass  # Already unlinked or not in scene collection
