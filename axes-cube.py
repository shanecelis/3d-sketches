import bpy
import colorsys
import math
import sys
from pathlib import Path
from mathutils import Matrix, Vector
from typing import Optional, Tuple

# ============================================================
# CLI ARG PARSING (Blender passes its own args; ours come after --)
# Usage:
#   blender --background --factory-startup --python axes-cube.py -- /path/to/axes-cube.glb
# ============================================================
def get_cli_output_path(default_name: str = "axes-cube.glb") -> str:
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
Path(out_path).parent.mkdir(parents=True, exist_ok=True)

# ============================================================
# CLEAN SCENE
# ============================================================
bpy.ops.wm.read_factory_settings(use_empty=True)
if bpy.context.mode != "OBJECT":
    bpy.ops.object.mode_set(mode="OBJECT")

# ============================================================
# PARAMETERS (tweak as needed)
# ============================================================
cube_size = 1.0

# Small visual separation between clickable parts (keep tiny; large values look "exploded").
gap = 0.002

# Face "panel" thickness along its normal.
face_thickness = 0.06

# Border edge beam thickness (square cross-section).
border_thickness = 0.08

# Corner cube size.
corner_size = 0.16

# Depth ordering to avoid z-fighting when parts get very close:
# corners sit slightly above borders, borders slightly above faces.
recess_face_length = 0.05
border_outset = 0.002
corner_outset = 0.004

# ============================================================
# HELPERS
# ============================================================
def ensure_collection(name: str, parent: Optional[bpy.types.Collection] = None) -> bpy.types.Collection:
    col = bpy.data.collections.get(name)
    if col is None:
        col = bpy.data.collections.new(name)
    if parent is None:
        # Blender 2.93: __contains__ expects a collection name (string), not a Collection object.
        if bpy.context.scene.collection.children.get(col.name) is None:
            bpy.context.scene.collection.children.link(col)
    else:
        if parent.children.get(col.name) is None:
            parent.children.link(col)
    return col


def make_mat(name: str, rgba: Tuple[float, float, float, float]) -> bpy.types.Material:
    mat = bpy.data.materials.get(name)
    if mat is None:
        mat = bpy.data.materials.new(name=name)
        mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf:
        bsdf.inputs["Base Color"].default_value = rgba
        bsdf.inputs["Roughness"].default_value = 0.6
        bsdf.inputs["Metallic"].default_value = 0.0
    return mat


def create_box(
    name: str,
    size_xyz: Tuple[float, float, float],
    location=(0.0, 0.0, 0.0),
    rotation_euler=(0.0, 0.0, 0.0),
):
    # Blender's "default cube" is 2x2x2. Our scaling below assumes that base size.
    bpy.ops.mesh.primitive_cube_add(size=2.0, location=location, rotation=rotation_euler)
    obj = bpy.context.object
    obj.name = name
    obj.scale = (size_xyz[0] / 2.0, size_xyz[1] / 2.0, size_xyz[2] / 2.0)
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
    return obj


def move_to_collection(obj: bpy.types.Object, col: bpy.types.Collection):
    for c in list(obj.users_collection):
        c.objects.unlink(obj)
    col.objects.link(obj)


def set_click_meta(obj: bpy.types.Object, part: str, label: str):
    # Makes it easy to distinguish parts in-engine (glTF exports these as extras).
    obj["click_part"] = part  # "face" | "border" | "corner"
    obj["click_label"] = label


def darken_rgba(rgba: Tuple[float, float, float, float], factor: float = 0.6) -> Tuple[float, float, float, float]:
    """Darken an RGBA color by scaling HSV value (brightness). factor in [0..1]."""
    r, g, b, a = rgba
    h, s, v = colorsys.rgb_to_hsv(r, g, b)
    v *= factor
    r2, g2, b2 = colorsys.hsv_to_rgb(h, s, v)
    return (r2, g2, b2, a)


# ============================================================
# COLLECTIONS
# ============================================================
root_col = ensure_collection("AxesCube")
faces_col = ensure_collection("Faces", root_col)
borders_col = ensure_collection("Borders", root_col)
corners_col = ensure_collection("Corners", root_col)

# ============================================================
# MATERIALS (distinct, so picking can also use material if desired)
# ============================================================
mat_face = make_mat("Mat_Face", (0.85, 0.85, 0.85, 1.0))
mat_border = make_mat("Mat_Border", (0.15, 0.15, 0.15, 1.0))
mat_corner = make_mat("Mat_Corner", (0.35, 0.35, 0.35, 1.0))
blue = (0.30, 0.70, 1.00, 1.0)
red = (1.00, 0.35, 0.35, 1.0)
green = (0.35, 1.00, 0.50, 1.0)
darken_factor = 0.3
dark_blue = darken_rgba(blue, darken_factor)
dark_red = darken_rgba(red, darken_factor)
dark_green = darken_rgba(green, darken_factor)
# Optional per-face tints (handy when debugging orientation)
face_tints = {
    "front": make_mat("Mat_Front", blue),  # +Y
    "right": make_mat("Mat_Right", red),  # +X
    "top": make_mat("Mat_Top", green),  # +Z
    "back": make_mat("Mat_Back", dark_blue),  # -Y
    "left": make_mat("Mat_Left", dark_red),  # -X
    "bottom": make_mat("Mat_Bottom", dark_green),  # -Z
}

# ============================================================
# BUILD PARTS
# Coordinate convention (Blender):
# - Front: +Y
# - Right: +X
# - Top: +Z
# - Back: -Y
# - Left: -X
# - Bottom: -Z
# ============================================================
half = cube_size / 2.0

# Face panels: inset so they don't fight with borders.
# (Corners/borders are kept distinct mostly by depth ordering + materials, not large gaps.)
inset = border_thickness + gap
panel_size = max(0.01, cube_size - 2.0 * inset)

faces = [
    ("front", Vector((0.0, 1.0, 0.0))),
    ("right", Vector((1.0, 0.0, 0.0))),
    ("top", Vector((0.0, 0.0, 1.0))),
    ("back", Vector((0.0, -1.0, 0.0))),
    ("left", Vector((-1.0, 0.0, 0.0))),
    ("bottom", Vector((0.0, 0.0, -1.0))),
]

for label, n in faces:
    center = n * (half - face_thickness / 2.0 - recess_face_length)

    # Default cube is aligned; rotate panel so its thickness axis matches face normal.
    if n == Vector((0.0, 1.0, 0.0)) or n == Vector((0.0, -1.0, 0.0)):
        size = (panel_size, face_thickness, panel_size)
        rot = (0.0, 0.0, 0.0)
    elif n == Vector((1.0, 0.0, 0.0)) or n == Vector((-1.0, 0.0, 0.0)):
        size = (face_thickness, panel_size, panel_size)
        rot = (0.0, 0.0, 0.0)
    else:
        size = (panel_size, panel_size, face_thickness)
        rot = (0.0, 0.0, 0.0)

    obj = create_box(f"Face_{label.capitalize()}", size, location=center, rotation_euler=rot)
    move_to_collection(obj, faces_col)
    obj.data.materials.append(face_tints.get(label, mat_face))
    set_click_meta(obj, "face", label)


# Border edges: 12 beams along cube edges, separated from corners.
edge_len = max(0.01, cube_size - 2.0 * corner_size - 2.0 * gap)
fixed = half - border_thickness / 2.0

def add_edge(name: str, size_xyz: tuple[float, float, float], loc: Vector):
    obj = create_box(name, size_xyz, location=loc, rotation_euler=(0.0, 0.0, 0.0))
    move_to_collection(obj, borders_col)
    obj.data.materials.append(mat_border)
    set_click_meta(obj, "border", name)


# Edges parallel to X (vary x; fixed y,z)
for sy in (-1.0, 1.0):
    for sz in (-1.0, 1.0):
        loc = Vector((0.0, sy * (fixed + border_outset), sz * (fixed + border_outset)))
        add_edge(f"Border_X_y{int(sy)}_z{int(sz)}", (edge_len, border_thickness, border_thickness), loc)

# Edges parallel to Y (vary y; fixed x,z)
for sx in (-1.0, 1.0):
    for sz in (-1.0, 1.0):
        loc = Vector((sx * (fixed + border_outset), 0.0, sz * (fixed + border_outset)))
        add_edge(f"Border_Y_x{int(sx)}_z{int(sz)}", (border_thickness, edge_len, border_thickness), loc)

# Edges parallel to Z (vary z; fixed x,y)
for sx in (-1.0, 1.0):
    for sy in (-1.0, 1.0):
        loc = Vector((sx * (fixed + border_outset), sy * (fixed + border_outset), 0.0))
        add_edge(f"Border_Z_x{int(sx)}_y{int(sy)}", (border_thickness, border_thickness, edge_len), loc)


# Corners: 8 cubes at vertices.
corner_fixed = half - corner_size / 2.0
for sx in (-1.0, 1.0):
    for sy in (-1.0, 1.0):
        for sz in (-1.0, 1.0):
            loc = Vector(
                (
                    sx * (corner_fixed + corner_outset),
                    sy * (corner_fixed + corner_outset),
                    sz * (corner_fixed + corner_outset),
                )
            )
            obj = create_box(f"Corner_x{int(sx)}_y{int(sy)}_z{int(sz)}", (corner_size, corner_size, corner_size), location=loc)
            move_to_collection(obj, corners_col)
            obj.data.materials.append(mat_corner)
            set_click_meta(obj, "corner", obj.name)


# Parent everything to an empty so you can move/scale as one.
axes_root = bpy.data.objects.new("AxesCubeRoot", None)
move_to_collection(axes_root, root_col)

for col in (faces_col, borders_col, corners_col):
    for obj in list(col.objects):
        if obj is axes_root:
            continue
        obj.parent = axes_root

# ============================================================
# EXPORT GLB
# ============================================================
bpy.ops.object.select_all(action="DESELECT")
axes_root.select_set(True)
for col in (faces_col, borders_col, corners_col):
    for obj in col.objects:
        obj.select_set(True)

bpy.context.view_layer.objects.active = axes_root

bpy.ops.export_scene.gltf(
    filepath=out_path,
    export_format="GLB",
    use_selection=True,
    export_yup=True,
    export_apply=False,
    export_extras=True,
    export_skins=False,
    export_animations=False,
)

print(f"Exported: {out_path}")

