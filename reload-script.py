#!/usr/bin/env python3
# Stick this in Blender's scripting window so you can easily reload/reimport
# from a file.

import bpy
import os

GLB_PATH = "/Users/shane/Projects/rocket-wave/assets/assemble-rocket.glb"
COLLECTION_NAME = "GLB_RELOAD"

# Delete previous import
if COLLECTION_NAME in bpy.data.collections:
    for obj in bpy.data.collections[COLLECTION_NAME].objects:
        bpy.data.objects.remove(obj, do_unlink=True)

# Import
bpy.ops.import_scene.gltf(filepath=GLB_PATH)

# Move imported objects into a known collection
col = bpy.data.collections.get(COLLECTION_NAME)
if not col:
    col = bpy.data.collections.new(COLLECTION_NAME)
    bpy.context.scene.collection.children.link(col)

for obj in bpy.context.selected_objects:
    col.objects.link(obj)
    bpy.context.scene.collection.objects.unlink(obj)
