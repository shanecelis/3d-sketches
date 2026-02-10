import bpy
import math
import sys
from pathlib import Path
from mathutils import Vector

# ============================================================
# CLI ARG PARSING (Blender passes its own args; ours come after --)
# Usage:
#   blender --background --factory-startup --python crazyflie-rig.py -- /path/to/out.glb [/path/to/in.glb]
# Defaults:
#   in:  <script_dir>/crazyflie.glb
#   out: <script_dir>/crazyflie-rigged.glb
# ============================================================


def get_user_args():
    argv = sys.argv
    if "--" in argv:
        return argv[argv.index("--") + 1 :]
    return []


def main():
    user_args = get_user_args()

    script_dir = Path(__file__).resolve().parent
    out_path = Path(user_args[0]).expanduser().resolve() if len(user_args) >= 1 else (script_dir / "crazyflie-rigged.glb")
    in_path = Path(user_args[1]).expanduser().resolve() if len(user_args) >= 2 else (script_dir / "crazyflie.glb")

    if not in_path.exists():
        raise RuntimeError(f"Missing input GLB: {in_path}")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Clean scene
    bpy.ops.wm.read_factory_settings(use_empty=True)
    if bpy.context.mode != "OBJECT":
        bpy.ops.object.mode_set(mode="OBJECT")

    # Import GLB
    before = set(bpy.data.objects)
    bpy.ops.import_scene.gltf(filepath=str(in_path))
    imported = [o for o in bpy.data.objects if o not in before]
    mesh_objs = [o for o in imported if o.type == "MESH"]
    if not mesh_objs:
        raise RuntimeError(f"No mesh objects imported from {in_path}")

    # Identify a stable root for transforms. We will make the armature the new root and parent all
    # imported objects under it while keeping world transforms. This preserves any scaling that was
    # previously inherited from glTF root nodes.
    imported_objs = [o for o in imported if o.name != ""]

    # ------------------------------------------------------------
    # Recenter the imported model so its combined mesh bounds are centered at the origin.
    # This changes the exported model "anchor" (origin) in-engine.
    # ------------------------------------------------------------
    def world_bounds(objs):
        pts = []
        for obj in objs:
            if obj.type != "MESH":
                continue
            for c in obj.bound_box:
                pts.append(obj.matrix_world @ Vector(c))
        if not pts:
            return None
        vmin = Vector((min(p.x for p in pts), min(p.y for p in pts), min(p.z for p in pts)))
        vmax = Vector((max(p.x for p in pts), max(p.y for p in pts), max(p.z for p in pts)))
        return vmin, vmax

    bounds = world_bounds(mesh_objs)
    if bounds is not None:
        vmin, vmax = bounds
        center = (vmin + vmax) * 0.5

        imported_set = set(imported_objs)
        roots = [o for o in imported_objs if (o.parent is None or o.parent not in imported_set)]
        for r in roots:
            mw = r.matrix_world.copy()
            mw.translation = mw.translation - center
            r.matrix_world = mw

        print(f"[crazyflie-rig] Recentering: moved roots by {-center}")
    else:
        print("[crazyflie-rig] Recentering: skipped (no mesh bounds found)")

    # Find propellers by name prefix
    props = sorted([o for o in mesh_objs if o.name.lower().startswith("propeller")], key=lambda o: o.name.lower())
    if not props:
        raise RuntimeError("No propeller meshes found (expected names starting with 'propeller').")

    print(f"[crazyflie-rig] Imported meshes: {len(mesh_objs)}; propellers: {len(props)}")

    # Create armature (start in Object mode so we can apply transforms safely)
    bpy.ops.object.armature_add(enter_editmode=False, location=(0.0, 0.0, 0.0))
    rig = bpy.context.object
    rig.name = "CrazyflieRig"
    arm = rig.data
    arm.name = "CrazyflieRigData"

    bpy.ops.object.mode_set(mode="EDIT")

    # Root bone (created by armature_add)
    root = arm.edit_bones[0]
    root.name = "Root"
    root.head = (0.0, 0.0, 0.0)
    root.tail = (0.0, 0.0, 0.2)

    # One bone per propeller: spin around prop local +Y
    prop_bone_names = []
    rig_inv = rig.matrix_world.inverted()
    for i, prop in enumerate(props):
        bone_name = f"Root.Propeller_{i}"
        head_world = prop.matrix_world.to_translation()
        head = rig_inv @ head_world
        axis_world = (prop.matrix_world.to_3x3() @ Vector((0.0, 1.0, 0.0))).normalized()
        axis = (rig_inv.to_3x3() @ axis_world).normalized()

        # Small tail length based on propeller bounds
        bb = [Vector(c) for c in prop.bound_box]
        extent = max((b - bb[0]).length for b in bb) if bb else 0.05
        length = max(0.03, min(0.15, extent * 0.35))
        tail = head + axis * length

        b = arm.edit_bones.new(bone_name)
        b.head = head
        b.tail = tail
        b.parent = root
        b.use_connect = False
        prop_bone_names.append(bone_name)

    bpy.ops.object.mode_set(mode="OBJECT")

    # Make the armature the new root: parent all imported objects to the rig while keeping transforms.
    # This prevents loss of inherited glTF scaling when exporting.
    for obj in imported_objs:
        if obj == rig:
            continue
        bpy.ops.object.select_all(action="DESELECT")
        obj.select_set(True)
        rig.select_set(True)
        bpy.context.view_layer.objects.active = rig
        bpy.ops.object.parent_set(type="OBJECT", keep_transform=True)

    # Rig propellers using armature deform (skinning) so glTF/Bevy handles the joint transforms.
    # Each propeller is 100% weighted to its own bone.
    for prop, bone_name in zip(props, prop_bone_names):
        # Ensure an Armature modifier exists and points to our rig.
        mod = prop.modifiers.get("Armature") or prop.modifiers.new(name="Armature", type="ARMATURE")
        mod.object = rig

        # Create vertex group matching bone name, weight all vertices to 1.0.
        vg = prop.vertex_groups.get(bone_name) or prop.vertex_groups.new(name=bone_name)
        idxs = [v.index for v in prop.data.vertices]
        vg.add(idxs, 1.0, "REPLACE")

    # Automated sanity check: rotating a prop bone should move propeller vertices (skinning).
    bpy.ops.object.select_all(action="DESELECT")
    bpy.context.view_layer.objects.active = rig
    bpy.ops.object.mode_set(mode="POSE")
    depsgraph = bpy.context.evaluated_depsgraph_get()

    for prop, bone_name in zip(props, prop_bone_names):
        pb = rig.pose.bones[bone_name]
        # Reset pose
        pb.rotation_mode = "XYZ"
        pb.rotation_euler = (0.0, 0.0, 0.0)
        depsgraph.update()
        prop_eval = prop.evaluated_get(depsgraph)
        mesh_eval = prop_eval.to_mesh()
        if not mesh_eval.vertices:
            prop_eval.to_mesh_clear()
            continue
        v0_before = prop_eval.matrix_world @ mesh_eval.vertices[0].co.copy()
        prop_eval.to_mesh_clear()

        # Rotate around local +Y
        pb.rotation_euler = (0.0, math.radians(15.0), 0.0)
        depsgraph.update()
        prop_eval = prop.evaluated_get(depsgraph)
        mesh_eval = prop_eval.to_mesh()
        v0_after = prop_eval.matrix_world @ mesh_eval.vertices[0].co.copy()
        prop_eval.to_mesh_clear()

        if (v0_after - v0_before).length < 1e-7:
            raise RuntimeError(f"Sanity check failed: prop '{prop.name}' did not deform/move when bone rotated.")

        # Reset
        pb.rotation_euler = (0.0, 0.0, 0.0)

    bpy.ops.object.mode_set(mode="OBJECT")
    print("[crazyflie-rig] Sanity check passed: prop meshes respond to bone rotations.")

    # Export GLB (select imported objects + rig)
    bpy.ops.object.select_all(action="DESELECT")
    for obj in mesh_objs:
        obj.select_set(True)
    rig.select_set(True)
    bpy.context.view_layer.objects.active = rig

    bpy.ops.export_scene.gltf(
        filepath=str(out_path),
        export_format="GLB",
        use_selection=True,
        export_yup=True,
        export_apply=False,
        export_skins=True,
        export_animations=False,
    )

    print(f"[crazyflie-rig] Exported: {out_path}")


if __name__ == "__main__":
    main()

