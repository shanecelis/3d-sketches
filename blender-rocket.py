import bpy
import math
import sys
from pathlib import Path
from mathutils import Vector, Matrix

# ============================================================
# CLI ARG PARSING (Blender passes its own args; ours come after --)
# Usage:
#   blender --background --factory-startup --python make_rocket_glb.py -- /path/to/rocket.glb
# ============================================================
def get_cli_output_path(default_name="rocket.glb") -> str:
    argv = sys.argv
    if "--" in argv:
        user_args = argv[argv.index("--") + 1 :]
    else:
        user_args = []

    if len(user_args) >= 1 and user_args[0].strip():
        return user_args[0]
    return str(Path.cwd() / default_name)

out_path = get_cli_output_path()
out_path = str(Path(out_path).expanduser().resolve())

# Ensure output directory exists
Path(out_path).parent.mkdir(parents=True, exist_ok=True)

# ============================================================
# CLEAN SCENE
# ============================================================
bpy.ops.wm.read_factory_settings(use_empty=True)

if bpy.context.mode != "OBJECT":
    bpy.ops.object.mode_set(mode="OBJECT")

scene = bpy.context.scene
scene.unit_settings.system = "METRIC"

# ============================================================
# PARAMETERS
# ============================================================
body_radius = 0.25
body_height = 2.0

cone_height = 0.6
cone_radius = body_radius

fin_height = 0.45
fin_length = 0.60
fin_thickness = 0.08

fin_count = 4
fin_z = -body_height / 2 + fin_height / 2 + 0.05

# Push fins outward so they are visible
fin_radial_offset = body_radius + fin_thickness * 0.6 + 0.01

# ============================================================
# HELPERS
# ============================================================
def weight_all_verts(obj, vg_name, weight=1.0):
    """Assign all vertices of obj to a vertex group."""
    vg = obj.vertex_groups.get(vg_name) or obj.vertex_groups.new(name=vg_name)
    idxs = [v.index for v in obj.data.vertices]
    vg.add(idxs, weight, "REPLACE")
    return vg

# ============================================================
# CREATE BODY
# ============================================================
bpy.ops.mesh.primitive_cylinder_add(
    radius=body_radius,
    depth=body_height,
    vertices=48,
    location=(0, 0, 0),
)
body = bpy.context.object
body.name = "Body"
weight_all_verts(body, "Root", 1.0)

# ============================================================
# CREATE NOSE CONE
# ============================================================
bpy.ops.mesh.primitive_cone_add(
    radius1=cone_radius,
    radius2=0.0,
    depth=cone_height,
    vertices=48,
    location=(0, 0, body_height / 2 + cone_height / 2),
)
cone = bpy.context.object
cone.name = "Nose"
weight_all_verts(cone, "Root", 1.0)

# ============================================================
# CREATE FINS (distributed around rocket)
# Each fin is fully weighted to Fin_i
# ============================================================
fins = []
for i in range(fin_count):
    angle = (2 * math.pi / fin_count) * i

    bpy.ops.mesh.primitive_cube_add(size=1.0)
    fin = bpy.context.object
    fin.name = f"FinObj_{i}"

    # Scale fin: X thickness, Y outward length, Z height
    fin.scale = (fin_height / 2, fin_thickness / 2, fin_length / 2)
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)

    outward_dir = Vector((math.cos(angle), math.sin(angle), 0.0))

    fin.location = outward_dir * fin_radial_offset
    fin.location.z = fin_z

    # Rotate so fin "points" outward
    fin.rotation_euler = (0.0, 0.0, angle)
    bpy.ops.object.transform_apply(location=False, rotation=True, scale=False)

    weight_all_verts(fin, f"Fin_{i}", 1.0)
    fins.append(fin)

# ============================================================
# JOIN INTO ONE MESH
# ============================================================
for obj in [body, cone] + fins:
    obj.select_set(True)
bpy.context.view_layer.objects.active = body
bpy.ops.object.join()

rocket = bpy.context.object
rocket.name = "Rocket"
bpy.ops.object.shade_smooth()

# ============================================================
# CREATE ARMATURE / RIG
# ============================================================
bpy.ops.object.armature_add(enter_editmode=True, location=(0, 0, 0))
rig = bpy.context.object
rig.name = "RocketRig"
arm = rig.data
arm.name = "RocketRigData"

# Root bone
root_bone = arm.edit_bones[0]
root_bone.name = "Root"
root_bone.head = (0, 0, -body_height / 2)
root_bone.tail = (0, 0, body_height / 2 + cone_height)

# Fin bones: hinge at (near) fin base, tail outward
for i in range(fin_count):
    angle = (2 * math.pi / fin_count) * i

    # Put hinge near body surface; bias slightly inward along outward axis
    hinge = Vector((fin_radial_offset - fin_length * 0.25, 0.0, fin_z))
    hinge = Matrix.Rotation(angle, 4, "Z") @ hinge

    outward = Vector((1.0, 0.0, 0.0))
    outward = Matrix.Rotation(angle, 4, "Z") @ outward

    b = arm.edit_bones.new(f"Fin_{i}")
    b.head = hinge
    b.tail = hinge + outward * (fin_length * 0.6)
    b.parent = root_bone
    b.use_connect = False

bpy.ops.object.mode_set(mode="OBJECT")

# Bind mesh to armature (skinning)
mod = rocket.modifiers.new(name="Armature", type="ARMATURE")
mod.object = rig
mod.use_vertex_groups = True
rocket.parent = rig

# Make rotation editing predictable (not required for export, but nice)
bpy.context.view_layer.objects.active = rig
bpy.ops.object.mode_set(mode="POSE")
for pb in rig.pose.bones:
    pb.rotation_mode = "XYZ"
bpy.ops.object.mode_set(mode="OBJECT")

# ============================================================
# EXPORT GLB
# ============================================================
# Select only our objects for export
bpy.ops.object.select_all(action="DESELECT")
rocket.select_set(True)
rig.select_set(True)
bpy.context.view_layer.objects.active = rocket

# Export options:
# - export_skins=True ensures the rig/weights go out
# - export_animations=False (no clips created here; you can animate in your engine)
bpy.ops.export_scene.gltf(
    filepath=out_path,
    export_format="GLB",
    use_selection=True,
    export_yup=True,
    export_apply=False,
    export_skins=True,
    export_animations=False,
)

print(f"Exported: {out_path}")
