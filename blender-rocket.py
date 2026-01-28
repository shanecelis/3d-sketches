import bpy
import math
from mathutils import Vector, Matrix

# -----------------------------
# Clean scene
# -----------------------------
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

if bpy.context.mode != 'OBJECT':
    bpy.ops.object.mode_set(mode='OBJECT')

scene = bpy.context.scene
scene.unit_settings.system = 'METRIC'

# -----------------------------
# Parameters
# -----------------------------
body_radius = 0.25
body_height = 2.0

cone_height = 0.6
cone_radius = body_radius

fin_height = 0.45
fin_length = 0.60
fin_thickness = 0.08

fin_z = -body_height / 2 + fin_height / 2 + 0.05
fin_count = 4

# Push fins outward so they don't get buried in the fuselage
fin_radial_offset = body_radius + fin_thickness * 0.6 + 0.01


# -----------------------------
# Helpers
# -----------------------------
def ensure_object_mode():
    if bpy.context.mode != 'OBJECT':
        bpy.ops.object.mode_set(mode='OBJECT')


def weight_all_verts(obj, vg_name, weight=1.0):
    """Assign all vertices of obj to a vertex group."""
    vg = obj.vertex_groups.get(vg_name) or obj.vertex_groups.new(name=vg_name)
    idxs = [v.index for v in obj.data.vertices]
    vg.add(idxs, weight, 'REPLACE')
    return vg


# -----------------------------
# Create rocket body (cylinder)
# -----------------------------
bpy.ops.mesh.primitive_cylinder_add(
    radius=body_radius,
    depth=body_height,
    vertices=48,
    location=(0, 0, 0)
)
body = bpy.context.object
body.name = "Body"
weight_all_verts(body, "Root", 1.0)

# -----------------------------
# Create nose cone
# -----------------------------
bpy.ops.mesh.primitive_cone_add(
    radius1=cone_radius,
    radius2=0.0,
    depth=cone_height,
    vertices=48,
    location=(0, 0, body_height / 2 + cone_height / 2)
)
cone = bpy.context.object
cone.name = "Nose"
weight_all_verts(cone, "Root", 1.0)

# -----------------------------
# Create fins (distributed around rocket)
# Each fin is a thin box, fully weighted to its own Fin_i vertex group.
# -----------------------------
fins = []

for i in range(fin_count):
    angle = (2 * math.pi / fin_count) * i

    # Create fin cube
    bpy.ops.mesh.primitive_cube_add(size=1.0)
    fin = bpy.context.object
    fin.name = f"FinObj_{i}"

    # Scale into fin shape (X outward length, Y thickness, Z height)
    fin.scale = (fin_height / 2, fin_thickness / 2, fin_length / 2)
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)

    # Outward direction around Z axis
    outward_dir = Vector((math.cos(angle), math.sin(angle), 0.0))

    # Place fin around the body near bottom
    fin.location = outward_dir * fin_radial_offset
    fin.location.z = fin_z

    # Rotate fin so its "length" points outward
    fin.rotation_euler = (0.0, 0.0, angle)
    bpy.ops.object.transform_apply(location=False, rotation=True, scale=False)

    # Weight to its bone group
    vg_name = f"Fin_{i}"
    weight_all_verts(fin, vg_name, 1.0)

    fins.append(fin)

# -----------------------------
# Join into one mesh (preserve vertex groups)
# -----------------------------
ensure_object_mode()
for obj in [body, cone] + fins:
    obj.select_set(True)
bpy.context.view_layer.objects.active = body
bpy.ops.object.join()

rocket = bpy.context.object
rocket.name = "Rocket"
bpy.ops.object.shade_smooth()

# -----------------------------
# Create Armature / Rig
# -----------------------------
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

# Fin bones: heads at hinge points on fuselage surface, tails pointing outward
for i in range(fin_count):
    angle = (2 * math.pi / fin_count) * i

    # Hinge point at body surface, aligned with fin position
    hinge = Vector((fin_radial_offset - fin_length * 0.25, 0.0, fin_z))
    hinge = Matrix.Rotation(angle, 4, 'Z') @ hinge

    # Outward direction for bone tail
    outward = Vector((1.0, 0.0, 0.0))
    outward = Matrix.Rotation(angle, 4, 'Z') @ outward

    b = arm.edit_bones.new(f"Fin_{i}")
    b.head = hinge
    b.tail = hinge + outward * (fin_length * 0.6)
    b.parent = root_bone
    b.use_connect = False

bpy.ops.object.mode_set(mode='OBJECT')

# -----------------------------
# Bind mesh to armature (skinning)
# -----------------------------
mod = rocket.modifiers.new(name="Armature", type='ARMATURE')
mod.object = rig
mod.use_vertex_groups = True

rocket.parent = rig

# Optional: don't force bones on top of mesh (set True if you want them visible)
rig.show_in_front = False

# Optional: make pose rotation easier to edit
bpy.context.view_layer.objects.active = rig
bpy.ops.object.mode_set(mode='POSE')
for pb in rig.pose.bones:
    pb.rotation_mode = 'XYZ'
bpy.ops.object.mode_set(mode='OBJECT')

print("Created Rocket + RocketRig. Animate by rotating Fin_0..Fin_3 bones in Pose Mode.")
