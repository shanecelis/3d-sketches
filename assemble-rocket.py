import bpy
import bmesh
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

# Optional CLI overrides after output path:
#   blender ... --python assemble-rocket.py -- out.glb [Body.glb] [Fin.glb] [FinGhost.glb] [--add-ghost-fins]
# Flags can appear anywhere after `--` (order-insensitive).
argv = sys.argv
if "--" in argv:
    user_args = argv[argv.index("--") + 1 :]
else:
    user_args = []

script_dir = Path(__file__).resolve().parent

# Split args into positional args (paths) and flags.
pos_args = []
flags = set()
for a in user_args:
    if isinstance(a, str) and a.startswith("--"):
        flags.add(a)
    else:
        pos_args.append(a)

body_glb_path = Path(pos_args[1]).expanduser().resolve() if len(pos_args) >= 2 else (script_dir / "Body.glb")
fin_glb_path = Path(pos_args[2]).expanduser().resolve() if len(pos_args) >= 3 else (script_dir / "Fin.glb")
fin_ghost_glb_path = (
    Path(pos_args[3]).expanduser().resolve()
    if len(pos_args) >= 4
    else (script_dir / "assets" / "FinGhost.glb")
)

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
fin_count = 4

# If true, also import FinGhost.glb and add additional joints named FinGhost_0..FinGhost_3.
# The ghost fin geometry is included in the exported GLB and fully weighted to those bones.
#
# Option behavior:
# - Pass `--add-ghost-fins` to enable explicitly.
# - `--add-ghost-fin` is accepted as a backwards-compatible alias.
# - Passing a 4th positional arg (FinGhost.glb path) also enables.
use_fin_ghost = "--add-ghost-fins" in flags


# Body import fixup (degrees). This rotation is baked into the final exported mesh.
body_y_rotation_degrees = 0.0

# Fins: if your `Fin.glb` is already positioned correctly relative to the rocket origin,
# we can just duplicate it and spin copies around the rocket axis.
rocket_axis = "X"  # "X" | "Y" | "Z" (axis that goes through the rocket body)
place_fins_from_template = True
fin_spin_axis = "X"  # axis to rotate copies around (usually same as rocket_axis)

# Where along the fin height the bone pivot (hinge) should be.
# 0.0 = bottom of fin, 1.0 = top of fin.
fin_anchor_t = 0.65

# (Only used when place_fins_from_template=False)
# Where along the rocket height to place the fins, as an absolute offset above the body's bottom.
# 2 inches = 0.0508 meters.
INCH_TO_M = 0.0254
fin_axis_offset_m = 2.0 * INCH_TO_M

# Move the rocket so the base is at the origin, then apply that translation into the hierarchy.
# Set to the distance along rocket_axis to shift (e.g. 0.15 if the base is at -0.15).
# The offset is applied and baked so the root stays at (0,0,0) in the export.
children_offset_along_axis = 1.2

# ============================================================
# HELPERS
# ============================================================
def weight_all_verts(obj, vg_name, weight=1.0):
    """Assign all vertices of obj to a vertex group."""
    vg = obj.vertex_groups.get(vg_name) or obj.vertex_groups.new(name=vg_name)
    idxs = [v.index for v in obj.data.vertices]
    vg.add(idxs, weight, "REPLACE")
    return vg


def make_mat(name: str, rgba):
    """Create or fetch a Principled BSDF material with the given base color."""
    mat = bpy.data.materials.get(name)
    if mat is None:
        mat = bpy.data.materials.new(name=name)
        mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf:
        bsdf.inputs["Base Color"].default_value = rgba
        bsdf.inputs["Roughness"].default_value = 0.65
        bsdf.inputs["Metallic"].default_value = 0.0
    return mat



def make_vertex_color_mat(name: str, layer_name: str):
    """
    Material that displays vertex colors (by name) as Base Color.
    Useful so importing the exported GLB back into Blender immediately shows the paint.
    """
    mat = bpy.data.materials.get(name)
    if mat is None:
        mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True

    nt = mat.node_tree
    nt.nodes.clear()

    out = nt.nodes.new("ShaderNodeOutputMaterial")
    out.location = (300, 0)
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    bsdf.location = (80, 0)
    attr = nt.nodes.new("ShaderNodeAttribute")
    attr.location = (-160, 0)
    attr.attribute_name = layer_name

    nt.links.new(attr.outputs["Color"], bsdf.inputs["Base Color"])
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])
    return mat


def make_image_texture_mat(name: str, image: bpy.types.Image, uv_layer_name: str, extension: str = "CLIP"):
    """
    Material that displays an Image Texture using the given UV map (reliable Blender preview).

    - extension="CLIP": clamp outside 0..1 (good for atlases)
    - extension="REPEAT": wrap (good for cylindrical body UVs, avoids seam interpolation artifacts)
    """
    mat = bpy.data.materials.get(name)
    if mat is None:
        mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True

    nt = mat.node_tree
    nt.nodes.clear()

    out = nt.nodes.new("ShaderNodeOutputMaterial")
    out.location = (420, 0)
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    bsdf.location = (200, 0)
    tex = nt.nodes.new("ShaderNodeTexImage")
    tex.location = (-60, 0)
    tex.image = image
    tex.interpolation = "Closest"
    tex.extension = extension

    uv = nt.nodes.new("ShaderNodeUVMap")
    uv.location = (-260, -40)
    uv.uv_map = uv_layer_name

    nt.links.new(uv.outputs["UV"], tex.inputs["Vector"])
    nt.links.new(tex.outputs["Color"], bsdf.inputs["Base Color"])
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])
    return mat


def make_image_texture_mat_with_alpha(
    name: str,
    image: bpy.types.Image,
    uv_layer_name: str,
    alpha: float,
    extension: str = "CLIP",
):
    """
    Image-texture material that uses the image's RGB, but forces transparency via a constant alpha.
    (Useful for ghost fins; exported as glTF alphaMode=BLEND.)
    """
    mat = bpy.data.materials.get(name)
    if mat is None:
        mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True

    nt = mat.node_tree
    nt.nodes.clear()

    out = nt.nodes.new("ShaderNodeOutputMaterial")
    out.location = (520, 0)
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    bsdf.location = (260, 0)

    tex = nt.nodes.new("ShaderNodeTexImage")
    tex.location = (-40, 0)
    tex.image = image
    tex.interpolation = "Closest"
    tex.extension = extension

    uv = nt.nodes.new("ShaderNodeUVMap")
    uv.location = (-240, -40)
    uv.uv_map = uv_layer_name

    alpha_val = nt.nodes.new("ShaderNodeValue")
    alpha_val.location = (80, -220)
    alpha_val.outputs[0].default_value = float(alpha)

    nt.links.new(uv.outputs["UV"], tex.inputs["Vector"])
    nt.links.new(tex.outputs["Color"], bsdf.inputs["Base Color"])
    nt.links.new(alpha_val.outputs[0], bsdf.inputs["Alpha"])
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])

    mat.blend_method = "BLEND"
    if hasattr(mat, "shadow_method"):
        mat.shadow_method = "NONE"  # Removed in Blender 4.3+
    return mat


def assign_material(obj, mat):
    obj.data.materials.clear()
    obj.data.materials.append(mat)


def assign_fin_two_sided_materials(fin_obj, local_axis: Vector, mat_a, mat_b, mat_edge):
    """
    Assign different materials to the two "broad" fin side faces.
    We detect those faces by checking if the face normal is aligned with ±local_axis.

    Important: This should be called while the fin's mesh is still in its local orientation
    (i.e. before applying any object rotation), so the axis is unambiguous.
    """
    axis = local_axis.normalized()

    fin_obj.data.materials.clear()
    fin_obj.data.materials.append(mat_a)  # slot 0
    fin_obj.data.materials.append(mat_b)  # slot 1
    fin_obj.data.materials.append(mat_edge)  # slot 2

    if bpy.context.mode != "OBJECT":
        bpy.ops.object.mode_set(mode="OBJECT")

    bm = bmesh.new()
    bm.from_mesh(fin_obj.data)
    bm.faces.ensure_lookup_table()

    # Only treat faces with normals ~parallel to axis as "sides".
    # Everything else gets the edge material.
    align_threshold = 0.92
    for f in bm.faces:
        d = f.normal.dot(axis)
        if abs(d) >= align_threshold:
            f.material_index = 0 if d > 0.0 else 1
        else:
            f.material_index = 2

    bm.to_mesh(fin_obj.data)
    bm.free()

# ============================================================
# IMPORT HELPERS
# ============================================================
def import_glb_mesh(filepath: Path, join_name: str):
    """Import a GLB and return a single joined mesh object (parents cleared, transforms applied)."""
    filepath = Path(filepath).expanduser().resolve()
    if not filepath.exists():
        raise RuntimeError(f"Missing GLB: {filepath}")

    before = set(bpy.data.objects)
    bpy.ops.import_scene.gltf(filepath=str(filepath))
    imported = [o for o in bpy.data.objects if o not in before]
    meshes = [o for o in imported if o.type == "MESH"]
    if not meshes:
        raise RuntimeError(f"No mesh objects found in: {filepath}")

    # Clear parenting (keep world transforms), then bake transforms into mesh data.
    for obj in meshes:
        bpy.ops.object.select_all(action="DESELECT")
        obj.select_set(True)
        bpy.context.view_layer.objects.active = obj
        if obj.parent is not None:
            bpy.ops.object.parent_clear(type="CLEAR_KEEP_TRANSFORM")
        bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)

    # Join into one mesh object for simpler skinning/export.
    bpy.ops.object.select_all(action="DESELECT")
    for obj in meshes:
        obj.select_set(True)
    bpy.context.view_layer.objects.active = meshes[0]
    bpy.ops.object.join()
    joined = bpy.context.object
    joined.name = join_name
    return joined


def bake_object_rotation(obj, rotation_euler_xyz):
    """Set an object's rotation_euler and bake it into the mesh data."""
    if bpy.context.mode != "OBJECT":
        bpy.ops.object.mode_set(mode="OBJECT")

    bpy.ops.object.select_all(action="DESELECT")
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj

    obj.rotation_euler = rotation_euler_xyz
    bpy.ops.object.transform_apply(location=False, rotation=True, scale=False)


def bake_object_translation(obj, translation_xyz):
    """Set an object's location and bake it (into mesh data / children)."""
    if bpy.context.mode != "OBJECT":
        bpy.ops.object.mode_set(mode="OBJECT")

    bpy.ops.object.select_all(action="DESELECT")
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj

    obj.location = translation_xyz
    bpy.ops.object.transform_apply(location=True, rotation=False, scale=False)


def local_aabb(obj):
    """Return (min_vec, max_vec) of the object's local-space AABB (after transforms applied)."""
    pts = [Vector(c) for c in obj.bound_box]
    vmin = Vector((min(p.x for p in pts), min(p.y for p in pts), min(p.z for p in pts)))
    vmax = Vector((max(p.x for p in pts), max(p.y for p in pts), max(p.z for p in pts)))
    return vmin, vmax


def axis_index(axis: str) -> int:
    axis = axis.upper()
    if axis == "X":
        return 0
    if axis == "Y":
        return 1
    if axis == "Z":
        return 2
    raise ValueError(f"Invalid axis: {axis}")


def radial_distance_from_axis(v: Vector, axis: str) -> float:
    axis = axis.upper()
    if axis == "X":
        return math.sqrt(v.y * v.y + v.z * v.z)
    if axis == "Y":
        return math.sqrt(v.x * v.x + v.z * v.z)
    if axis == "Z":
        return math.sqrt(v.x * v.x + v.y * v.y)
    raise ValueError(f"Invalid axis: {axis}")


def estimate_hinge_point(mesh_obj, axis: str, t: float) -> Vector:
    """
    Estimate a good hinge point for the fin from its geometry:
    - choose vertices near the requested axis position (based on t along the fin extent)
    - then choose the ones closest to the rocket axis (minimum radial distance)
    """
    axis = axis.upper()
    ai = axis_index(axis)

    vmin, vmax = local_aabb(mesh_obj)
    axis_min = vmin[ai]
    axis_max = vmax[ai]
    axis_len = max(1e-6, axis_max - axis_min)
    target = axis_min + axis_len * float(t)
    eps = axis_len * 0.08

    verts = [v.co.copy() for v in mesh_obj.data.vertices]
    candidates = [v for v in verts if abs(v[ai] - target) <= eps]
    if not candidates:
        candidates = verts

    dists = [radial_distance_from_axis(v, axis) for v in candidates]
    min_d = min(dists)
    # take a small band near the minimum to stabilize average
    band = min_d + max(1e-6, 0.02 * axis_len)
    near = [v for v, d in zip(candidates, dists) if d <= band]
    if not near:
        near = candidates

    acc = Vector((0.0, 0.0, 0.0))
    for v in near:
        acc += v
    return acc / len(near)


def outward_dir_from_hinge(hinge: Vector, axis: str) -> Vector:
    axis = axis.upper()
    if axis == "X":
        v = Vector((0.0, hinge.y, hinge.z))
    elif axis == "Y":
        v = Vector((hinge.x, 0.0, hinge.z))
    else:  # "Z"
        v = Vector((hinge.x, hinge.y, 0.0))
    if v.length < 1e-6:
        return Vector((0.0, 1.0, 0.0))
    return v.normalized()


def hinge_pivot_quat(hinge_axis_world: Vector, outward_world: Vector):
    """
    Build a pivot orientation such that:
    - local +Y is the hinge axis (1-DOF rotation axis)
    - local +Z points roughly outward from the rocket body
    This makes it easy to animate by rotating about a single local axis.
    """
    y = hinge_axis_world.normalized()
    z = outward_world - y * outward_world.dot(y)
    if z.length < 1e-6:
        # Fallback if outward is accidentally parallel to hinge axis.
        z = Vector((0.0, 0.0, 1.0)) if abs(y.dot(Vector((0.0, 0.0, 1.0)))) < 0.9 else Vector((0.0, 1.0, 0.0))
        z = z - y * z.dot(y)
    z.normalize()
    x = y.cross(z).normalized()

    # Construct a rotation matrix with columns [x y z] in world space.
    m = Matrix.Identity(3)
    m.col[0] = x
    m.col[1] = y
    m.col[2] = z
    return m.to_quaternion()


def axis_radial_max(mesh_obj, axis: str) -> float:
    vals = [radial_distance_from_axis(v.co, axis) for v in mesh_obj.data.vertices]
    return max(vals) if vals else 0.0


def ensure_uv_layer(mesh: bpy.types.Mesh, uv_name: str):
    """Ensure a UV map exists (Blender join/export tends to keep only active object's UV maps)."""
    uv = mesh.uv_layers.get(uv_name)
    if uv is None:
        uv = mesh.uv_layers.new(name=uv_name)
    mesh.uv_layers.active = uv
    mesh.uv_layers.active_index = mesh.uv_layers.find(uv.name)
    return uv


def draw_disc(pixels: list, size: int, cx: float, cy: float, r_px: float, rgba):
    """
    Draw a filled circle into a flat RGBA float pixel buffer (row-major, size x size).
    Use for marker dots etc.; replace with other shapes as needed.
    """
    r2 = r_px * r_px
    x0 = max(0, int(cx - r_px))
    x1 = min(size - 1, int(cx + r_px))
    y0 = max(0, int(cy - r_px))
    y1 = min(size - 1, int(cy + r_px))
    for yy in range(y0, y1 + 1):
        dy = yy - cy
        for xx in range(x0, x1 + 1):
            dx = xx - cx
            if dx * dx + dy * dy <= r2:
                idx = (yy * size + xx) * 4
                pixels[idx : idx + 4] = list(rgba)


def draw_line(
    pixels: list,
    size: int,
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    rgba,
    width: int = 1,
):
    """
    Draw a line into a flat RGBA float pixel buffer (row-major, size x size).
    Uses Bresenham's line algorithm; coordinates are clipped to image bounds.
    width: line thickness in pixels (default 1); drawn by offsetting along the perpendicular.
    """
    if width <= 1:
        _draw_line_1px(pixels, size, x0, y0, x1, y1, rgba)
        return
    # Thick line: unit perpendicular offset, then draw multiple 1-pixel lines.
    dx = x1 - x0
    dy = y1 - y0
    length = math.sqrt(dx * dx + dy * dy)
    if length < 1e-9:
        _draw_line_1px(pixels, size, x0, y0, x1, y1, rgba)
        return
    px = -dy / length
    py = dx / length
    half = (width - 1) / 2.0
    for k in range(width):
        offset = (k - half) * 1.0
        _draw_line_1px(
            pixels, size,
            x0 + px * offset, y0 + py * offset,
            x1 + px * offset, y1 + py * offset,
            rgba,
        )


def _draw_line_1px(pixels: list, size: int, x0: float, y0: float, x1: float, y1: float, rgba):
    """Bresenham 1-pixel line; used by draw_line."""
    ix0, iy0 = int(round(x0)), int(round(y0))
    ix1, iy1 = int(round(x1)), int(round(y1))
    dx = abs(ix1 - ix0)
    dy = -abs(iy1 - iy0)
    sx = 1 if ix0 < ix1 else -1
    sy = 1 if iy0 < iy1 else -1
    err = dx + dy
    rgba_list = list(rgba)
    while True:
        if 0 <= ix0 < size and 0 <= iy0 < size:
            idx = (iy0 * size + ix0) * 4
            pixels[idx : idx + 4] = rgba_list
        if ix0 == ix1 and iy0 == iy1:
            break
        e2 = 2 * err
        if e2 >= dy:
            err += dy
            ix0 += sx
        if e2 <= dx:
            err += dx
            iy0 += sy


def draw_line_marker(pixels: list, size: int, cx: float, cy: float, r_px: float, rgba):
    """
    Draw a line-based marker with the same signature as draw_disc (pixels, size, cx, cy, r_px, rgba).
    Draws a plus/cross centered at (cx, cy) with arms of length r_px.
    """
    # draw_line(pixels, size, cx - r_px, cy, cx + r_px, cy, rgba, 3)
    draw_line(pixels, size, cx, cy - r_px, cx, cy + r_px, rgba, 5)


# Default marker used by body paint and fin atlas; set to draw_line_marker to use cross markers.
draw_marker = draw_disc
# draw_marker = draw_line_marker


def create_body_paint_image_with_cap(
    name: str,
    size: int,
    side_band_v0: float,
    side_band_v1: float,
    cap_v0: float,
    cap_v1: float,
    sector_colors,
):
    """
    Create an image where:
    - everything is black
    - side band (v in [side_band_v0..side_band_v1]) is split into 4 U sectors
    - bottom cap region (v in [cap_v0..cap_v1]) is a 4-quadrant disc
    """
    img = bpy.data.images.get(name)
    if img is None:
        img = bpy.data.images.new(name=name, width=size, height=size, alpha=True, float_buffer=True)
    else:
        img.scale(size, size)

    pixels = [0.0] * (size * size * 4)

    def set_px(x, y, rgba):
        idx = (y * size + x) * 4
        pixels[idx : idx + 4] = list(rgba)

    def draw_dots(cx, cy, dot_count: int, spacing_px: int, r_outer: int):
        """Draw 1/2/3 dots as white ring + black center for visibility."""
        if dot_count <= 0:
            return
        r_inner = max(1, int(r_outer * 0.55))
        start_x = cx - (dot_count - 1) * spacing_px // 2
        for k in range(dot_count):
            x = start_x + k * spacing_px
            draw_marker(pixels, size, x, cy, r_inner, (0.0, 0.0, 0.0, 1.0))

    black = (0.0, 0.0, 0.0, 1.0)

    # Fill all black.
    for y in range(size):
        for x in range(size):
            set_px(x, y, black)

    def y_from_v(v: float) -> int:
        # Our pixels are written with y=0 at bottom of the image.
        return max(0, min(size, int(round(size * max(0.0, min(1.0, v))))))

    # Side band: 4 sectors by U.
    y0 = y_from_v(side_band_v0)
    y1 = y_from_v(side_band_v1)
    if y1 < y0:
        y0, y1 = y1, y0
    if y1 == y0:
        y1 = min(size, y0 + 1)

    for y in range(y0, y1):
        for x in range(size):
            sector = min(3, int((x / size) * 4))
            set_px(x, y, sector_colors[sector])

    # Bottom cap: quadrant disc in [cap_v0..cap_v1].
    # IMPORTANT: represent the cap as a true circle in a SQUARE sub-rect of the texture,
    # otherwise it will be stretched into an ellipse (since the cap is only a V slice).
    cy0 = y_from_v(cap_v0)
    cy1 = y_from_v(cap_v1)
    if cy1 < cy0:
        cy0, cy1 = cy1, cy0
    if cy1 == cy0:
        cy1 = min(size, cy0 + 1)

    cap_h = max(1, cy1 - cy0)
    # Square region for the cap disc: width == height == cap_h, centered in U.
    cap_x0 = (size - cap_h) // 2
    cap_x1 = min(size, cap_x0 + cap_h)
    for y in range(cy0, cy1):
        # local v in [0..1]
        lv = (y - cy0 + 0.5) / cap_h
        for x in range(cap_x0, cap_x1):
            u_local = (x - cap_x0 + 0.5) / cap_h
            yy = (u_local - 0.5) * 2.0
            zz = (lv - 0.5) * 2.0
            r = math.sqrt(yy * yy + zz * zz)
            if r > 1.0:
                continue  # keep black outside the disc
            ang = math.atan2(zz, yy)
            if ang < 0:
                ang += 2.0 * math.pi
            sector = int(ang / (math.pi / 2.0)) % 4
            set_px(x, y, sector_colors[sector])

    # Add fin index markers on the bottom cap, aligned to the quadrant axes (in texture space):
    # - Fin_1: along -X axis (left from center)
    # - Fin_2: along -Y axis (down from center)
    # - Fin_3: along +X axis (right from center)
    #
    # Dots are colinear on those axes (not a horizontal row).
    cap_r = 0.5 * cap_h
    r_start = 0.90 * cap_r
    r_step = 0.18 * cap_r
    cap_dot_r = max(2, int(cap_r * 0.12))

    cx0 = cap_x0 + int(round(cap_h * 0.5))
    cy_center = cy0 + int(round(cap_h * 0.5))

    def draw_dots_along(cx, cy, dot_count: int, dir_x: int, dir_y: int):
        for k in range(dot_count):
            if k >= 2:
                k += 1
            rr = r_start - k * r_step
            x = int(round(cx + dir_x * rr))
            y = int(round(cy + dir_y * rr))
            draw_marker(pixels, size, x, y, max(1, int(cap_dot_r * 0.55)), (0.0, 0.0, 0.0, 1.0))

    # Fin_1: 1 dot on -X
    draw_dots_along(cx0, cy_center, 1, -1, 0)
    # Fin_2: 2 dots on -Y
    draw_dots_along(cx0, cy_center, 2, 0, -1)
    # Fin_3: 3 dots on +X
    draw_dots_along(cx0, cy_center, 3, 1, 0)

    img.pixels = pixels
    try:
        img.pack()
    except Exception:
        pass
    img.colorspace_settings.name = "sRGB"
    return img


def composite_logo_into_body_paint(
    body_img: bpy.types.Image,
    logo_img: bpy.types.Image,
    u_center: float,
    v_center: float,
    logo_height_v: float,
    logo_rotation_degrees: float,
):
    """
    Alpha-composite `logo_img` onto `body_img` in UV space.
    - U wraps (logo can cross seam)
    - V is clamped
    """
    size = int(body_img.size[0])
    if size <= 0 or body_img.size[1] != size:
        raise RuntimeError("Body paint image must be square")

    lw, lh = int(logo_img.size[0]), int(logo_img.size[1])
    if lw <= 0 or lh <= 0:
        raise RuntimeError("Logo image has invalid dimensions")

    # Read pixels (flat RGBA float arrays).
    dst = list(body_img.pixels)
    src = list(logo_img.pixels)

    # Desired logo size in pixels (preserve aspect ratio).
    logo_h_px = max(1, int(round(size * max(0.0, min(1.0, logo_height_v)))))
    logo_w_px = max(1, int(round(logo_h_px * (lw / float(lh)))))

    x_center = int(round((u_center % 1.0) * size))
    y_center = int(round(max(0.0, min(1.0, v_center)) * size))
    x0 = x_center - logo_w_px // 2
    y0 = y_center - logo_h_px // 2

    def get_src_rgba(sx: int, sy: int):
        sx = max(0, min(lw - 1, sx))
        sy = max(0, min(lh - 1, sy))
        idx = (sy * lw + sx) * 4
        return src[idx], src[idx + 1], src[idx + 2], src[idx + 3]

    def blend(dst_rgba, src_rgba):
        dr, dg, db, da = dst_rgba
        sr, sg, sb, sa = src_rgba
        out_a = sa + da * (1.0 - sa)
        if out_a <= 1e-6:
            return 0.0, 0.0, 0.0, 0.0
        out_r = (sr * sa + dr * da * (1.0 - sa)) / out_a
        out_g = (sg * sa + dg * da * (1.0 - sa)) / out_a
        out_b = (sb * sa + db * da * (1.0 - sa)) / out_a
        return out_r, out_g, out_b, out_a

    theta = math.radians(logo_rotation_degrees)
    c, s = math.cos(theta), math.sin(theta)

    for yy in range(logo_h_px):
        y = y0 + yy
        if y < 0 or y >= size:
            continue

        for xx in range(logo_w_px):
            x = x0 + xx
            xw = x % size  # wrap around U

            # Logo-local coords in [-0.5, 0.5], rotate around center.
            u = (xx + 0.5) / logo_w_px
            v = (yy + 0.5) / logo_h_px
            lx = u - 0.5
            ly = v - 0.5
            rx = lx * c - ly * s
            ry = lx * s + ly * c
            su = rx + 0.5
            sv = ry + 0.5
            sx = int(su * lw)
            sy = int(sv * lh)

            sr, sg, sb, sa = get_src_rgba(sx, sy)
            if sa <= 1e-4:
                continue

            didx = (y * size + xw) * 4
            dr, dg, db, da = dst[didx], dst[didx + 1], dst[didx + 2], dst[didx + 3]
            out = blend((dr, dg, db, da), (sr, sg, sb, sa))
            dst[didx : didx + 4] = out

    body_img.pixels = dst


def save_image_png(img: bpy.types.Image, png_path: Path):
    """
    Save a Blender image to PNG, and verify the file exists on disk.
    In some Blender/background configurations, Image.save() can silently fail to produce a file,
    so we fall back to save_render().
    """
    png_path = Path(png_path).expanduser().resolve()
    png_path.parent.mkdir(parents=True, exist_ok=True)

    img.filepath_raw = str(png_path)
    img.file_format = "PNG"
    try:
        img.save()
    except Exception as e:
        print(f"[assemble-rocket] WARNING: Image.save() failed for {png_path}: {e}")

    if not png_path.exists():
        try:
            img.save_render(filepath=str(png_path))
        except Exception as e:
            print(f"[assemble-rocket] WARNING: Image.save_render() failed for {png_path}: {e}")

    if not png_path.exists():
        print(f"[assemble-rocket] WARNING: PNG was not created on disk: {png_path}")


def assign_body_uvs_with_bottom_cap_quadrants(
    mesh_obj,
    axis: str,
    uv_layer_name: str,
    u_offset_degrees: float,
    side_v_max: float,
    cap_v0: float,
    cap_v1: float,
    cap_rotation_degrees: float = 0.0,
    cap_uv_zoom_out: float = 1.0,
):
    """
    Assign UVs for the body:
    - side faces: cylindrical projection (U=angle, V=axis) mapped into v in [0..side_v_max]
    - bottom cap (negative axis normal): planar projection into cap region [cap_v0..cap_v1]
    - top cap: constant UV (black)
    - cap_uv_zoom_out: scale cap UVs so rim maps further toward circle edge (>1 = zoom out, less pixelation)
    """
    axis = axis.upper()
    if bpy.context.mode != "OBJECT":
        bpy.ops.object.mode_set(mode="OBJECT")

    mesh = mesh_obj.data
    bm = bmesh.new()
    bm.from_mesh(mesh)
    bm.faces.ensure_lookup_table()
    uv_layer = bm.loops.layers.uv.get(uv_layer_name)
    if uv_layer is None:
        uv_layer = bm.loops.layers.uv.new(uv_layer_name)

    vmin, vmax = local_aabb(mesh_obj)
    if axis == "X":
        amin, amax = vmin.x, vmax.x
        # radius is in YZ plane
        rmax = max(math.sqrt(v.co.y * v.co.y + v.co.z * v.co.z) for v in bm.verts) if bm.verts else 1.0
        axis_vec = Vector((1.0, 0.0, 0.0))
    elif axis == "Y":
        amin, amax = vmin.y, vmax.y
        rmax = max(math.sqrt(v.co.x * v.co.x + v.co.z * v.co.z) for v in bm.verts) if bm.verts else 1.0
        axis_vec = Vector((0.0, 1.0, 0.0))
    else:
        amin, amax = vmin.z, vmax.z
        rmax = max(math.sqrt(v.co.x * v.co.x + v.co.y * v.co.y) for v in bm.verts) if bm.verts else 1.0
        axis_vec = Vector((0.0, 0.0, 1.0))

    alen = max(1e-6, amax - amin)
    rmax = max(1e-6, rmax)

    # Use the bottom cap's actual radial extent so the cap rim maps to the circle edge
    # (otherwise rmax from the whole mesh can be larger than the cap, causing "zoomed in" UVs).
    cap_threshold = 0.9
    cap_rmax = 1e-6
    for f in bm.faces:
        if f.normal.dot(axis_vec) <= -cap_threshold:
            for v in f.verts:
                if axis == "X":
                    r = math.sqrt(v.co.y * v.co.y + v.co.z * v.co.z)
                elif axis == "Y":
                    r = math.sqrt(v.co.x * v.co.x + v.co.z * v.co.z)
                else:
                    r = math.sqrt(v.co.x * v.co.x + v.co.y * v.co.y)
                cap_rmax = max(cap_rmax, r)
    if cap_rmax <= 1e-6:
        cap_rmax = rmax

    two_pi = 2.0 * math.pi
    u_offset = (u_offset_degrees / 360.0) % 1.0
    cap_rot = Matrix.Rotation(math.radians(u_offset_degrees + cap_rotation_degrees), 4, axis)

    for f in bm.faces:
        nd = f.normal.dot(axis_vec)
        is_bottom_cap = nd <= -cap_threshold
        is_top_cap = nd >= cap_threshold

        for loop in f.loops:
            co = loop.vert.co

            if is_top_cap:
                loop[uv_layer].uv = (0.5, 0.5)  # black region
                continue

            if is_bottom_cap:
                # Rotate the radial plane so the quadrant pattern aligns with fins.
                co_r = cap_rot @ co
                if axis == "X":
                    yy, zz = co_r.y, co_r.z
                elif axis == "Y":
                    yy, zz = co_r.x, co_r.z
                else:
                    yy, zz = co_r.x, co_r.y

                # Map the cap disc into a SQUARE sub-rect of UV space:
                # - V span is [cap_v0..cap_v1] (height = cap_s)
                # - U span is centered with the same width cap_s (so circles stay circles)
                # cap_uv_zoom_out > 1 scales UVs outward so the rim reaches the circle edge (less pixelation).
                cap_s = max(1e-6, cap_v1 - cap_v0)
                cap_u0 = 0.5 - cap_s * 0.5

                u_local = 0.5 + (yy / (2.0 * cap_rmax)) * cap_uv_zoom_out
                v_local = 0.5 + (zz / (2.0 * cap_rmax)) * cap_uv_zoom_out
                u = cap_u0 + u_local * cap_s
                v = cap_v0 + v_local * cap_s
                loop[uv_layer].uv = (max(0.0, min(1.0, u)), max(0.0, min(1.0, v)))
                continue

            # Side faces: cylindrical.
            # We compute U in [0..1), then "unwrap" per-face so triangles crossing the seam
            # don't interpolate through the other sectors (the artifact you saw).
            if axis == "X":
                u0 = (math.atan2(co.z, co.y) / two_pi) + 0.5
                v0 = (co.x - amin) / alen
            elif axis == "Y":
                u0 = (math.atan2(co.x, co.z) / two_pi) + 0.5
                v0 = (co.y - amin) / alen
            else:  # "Z"
                u0 = (math.atan2(co.y, co.x) / two_pi) + 0.5
                v0 = (co.z - amin) / alen

            u_mod = (u0 + u_offset) % 1.0
            v = max(0.0, min(1.0, v0)) * side_v_max
            loop[uv_layer].uv = (u_mod, v)

        # Seam unwrap pass for side faces: if a face spans both ends of [0..1), shift low-U up by +1.
        if not (is_bottom_cap or is_top_cap):
            uvs = [loop[uv_layer].uv.copy() for loop in f.loops]
            us = [uv.x for uv in uvs]
            if us and (max(us) - min(us) > 0.5):
                for loop in f.loops:
                    uv = loop[uv_layer].uv
                    if uv.x < 0.5:
                        uv.x += 1.0
                        loop[uv_layer].uv = uv

    bm.to_mesh(mesh)
    bm.free()


def detect_fin_side_axis(mesh_obj, align_threshold: float = 0.80) -> Vector:
    """
    Pick which local axis (X/Y/Z) best matches the fin's broad side face normals.
    Returns a unit Vector in local space.
    """
    mesh = mesh_obj.data
    mesh.calc_normals()
    candidates = {
        "X": Vector((1.0, 0.0, 0.0)),
        "Y": Vector((0.0, 1.0, 0.0)),
        "Z": Vector((0.0, 0.0, 1.0)),
    }
    best = "Y"
    best_score = -1
    for k, axis in candidates.items():
        score = sum(1 for p in mesh.polygons if abs(p.normal.dot(axis)) >= align_threshold)
        if score > best_score:
            best = k
            best_score = score
    return candidates[best]


def detect_fin_thickness_axis_by_aabb(mesh_obj) -> Vector:
    """
    Detect fin thickness axis by choosing the smallest local AABB extent.
    This is usually the best proxy for "the two broad side faces are ±this axis".
    """
    vmin, vmax = local_aabb(mesh_obj)
    extents = {
        "X": vmax.x - vmin.x,
        "Y": vmax.y - vmin.y,
        "Z": vmax.z - vmin.z,
    }
    axis = min(extents.items(), key=lambda kv: kv[1])[0]
    return {"X": Vector((1.0, 0.0, 0.0)), "Y": Vector((0.0, 1.0, 0.0)), "Z": Vector((0.0, 0.0, 1.0))}[axis]

# ============================================================
# MATERIALS
# ============================================================
mat_body = make_mat("Mat_Body", (0.75, 0.75, 0.78, 1.0))
# Fin edge material (shared across fins).
mat_fin_edge = make_mat("Mat_Fin_Edge", (0.18, 0.18, 0.20, 1.0))

black = (0.0, 0.0, 0.0, 1.0)
red = (1.0, 0.37, 0.43, 1.0)
green = (0.29, 0.82, 0.44, 1.0)
yellow = (0.86, 0.93, 0.44, 1.0)
orange = (0.99, 0.71, 0.40, 1.0)
# Per-fin side colors.
# Customize these to whatever you want: 4 fins × (outer, inner).
fin_outer_colors = [
    # (1.00, 0.35, 0.35, 1.0),  # Fin 0 outer
    red,  # Fin 0 outer
    green,  # Fin 1 outer
    orange,  # Fin 2 outer
    yellow,  # Fin 3 outer
]

fin_inner_colors = [
    # (0.30, 0.70, 1.00, 1.0),  # Fin 0 inner
    green,  # Fin 0 inner
    orange,  # Fin 1 inner
    yellow,  # Fin 2 inner
    red,  # Fin 3 inner
]

# Body paint via texture (black body + colored band + colored bottom cap quadrants).
use_body_paint_texture = True
body_paint_image_name = "BodyPaint"
body_paint_size = 1024  # px
body_uv_layer = "UVMap"

# Logo (placed on the cylinder side region, i.e. before the bottom cap area).
use_body_logo = True
body_logo_path = Path(__file__).resolve().parent / "logo 10-8.png"
body_logo_u_center_degrees = 180.0 - 35.0  # 0° points at +Y when rocket_axis=="X"
body_logo_v_center_t = 0.55  # 0..1 along the side region only
body_logo_height_t = 0.18  # 0..1 of the side region height (keeps aspect ratio)
body_logo_rotation_degrees = 90.0

# Band is measured from the body's "bottom" along rocket_axis (min axis value).
body_band_start_m = 0.0
body_band_length_m = 0.0508  # 2 inches

# Rotate the colored sectors around the rocket axis (degrees) to align with fins.
body_sector_u_offset_degrees = 135.0

# Reserve top portion of the texture for the bottom-cap quadrant disc.
body_cap_v0 = 0.75
body_cap_v1 = 1.0
body_cap_rotation_degrees = 180.0
# Scale bottom-cap UVs so the rim maps to the circle edge (>1 = zoom out, reduces pixelation of markers).
body_cap_uv_zoom_out = 2.3

write_body_paint_png = True

# Paint fin sides via vertex colors (glTF COLOR_0) instead of using multiple materials.
# This is often nicer for engines: fewer primitives, but still distinct colors per side.
use_vertex_colors_for_fin_sides = False
fin_edge_color = (0.18, 0.18, 0.20, 1.0)
fin_vertex_color_layer = "Col"

# Alternative: paint fin sides via a generated texture atlas + UVs.
# This previews reliably in Blender after import.
use_fin_texture_atlas = True
fin_atlas_image_name = "FinAtlas"
fin_atlas_size = 512  # px
fin_uv_layer = "UVMap"

# Which local axis represents the fin "side" normal (the two broad faces are ±this axis).
# "AUTO" tries X/Y/Z and picks the one that best matches face normals.
fin_side_axis_mode = "AUTO"  # "AUTO" | "X" | "Y" | "Z"
fin_side_align_threshold = 0.65

# Write the generated fin atlas to disk so it can be inspected.
write_fin_atlas_png = True

# ============================================================
# LOAD BODY + FIN TEMPLATE
# ============================================================
body = import_glb_mesh(body_glb_path, "Body")
fin_template = import_glb_mesh(fin_glb_path, "FinTemplate")
fin_ghost_template = None

# Rotate the body to the desired orientation and bake it in.
bake_object_rotation(body, (0.0, math.radians(body_y_rotation_degrees), 0.0))

# If you want to preserve imported body materials, comment out the next line.
assign_material(body, mat_body)

# Ensure the active object's UV map exists before we join meshes later.
# (Join tends to preserve UV maps from the active object; we keep names consistent.)
ensure_uv_layer(body.data, body_uv_layer)

body_min, body_max = local_aabb(body)
ai = axis_index(rocket_axis)
body_axis_min = body_min[ai]
body_axis_max = body_max[ai]

# Paint the body (black everywhere, colored band + bottom cap quadrants).
if use_body_paint_texture:
    side_v_max = body_cap_v0
    body_axis_len = max(1e-6, body_axis_max - body_axis_min)
    band_v0 = (body_band_start_m / body_axis_len) * side_v_max
    band_v1 = ((body_band_start_m + body_band_length_m) / body_axis_len) * side_v_max

    assign_body_uvs_with_bottom_cap_quadrants(
        body,
        rocket_axis,
        body_uv_layer,
        u_offset_degrees=body_sector_u_offset_degrees,
        side_v_max=side_v_max,
        cap_v0=body_cap_v0,
        cap_v1=body_cap_v1,
        cap_rotation_degrees=body_cap_rotation_degrees,
        cap_uv_zoom_out=body_cap_uv_zoom_out,
    )

    body_img = create_body_paint_image_with_cap(
        body_paint_image_name,
        body_paint_size,
        side_band_v0=band_v0,
        side_band_v1=band_v1,
        cap_v0=body_cap_v0,
        cap_v1=body_cap_v1,
        sector_colors=fin_outer_colors,
    )

    # Composite the logo onto the side region (not the bottom cap).
    if use_body_logo and body_logo_path.exists():
        try:
            #print(f"[assemble-rocket] WARNING: loading rocket body image from {body_logo_path}")
            logo_img = bpy.data.images.load(filepath=str(body_logo_path), check_existing=True)
            # Place logo using side-region UVs.
            side_v_max = body_cap_v0
            u_center = (body_logo_u_center_degrees / 360.0 + 0.5) % 1.0
            v_center = max(0.0, min(1.0, body_logo_v_center_t)) * side_v_max
            logo_h_v = max(0.0, min(1.0, body_logo_height_t)) * side_v_max
            composite_logo_into_body_paint(
                body_img,
                logo_img,
                u_center=u_center,
                v_center=v_center,
                logo_height_v=logo_h_v,
                logo_rotation_degrees=body_logo_rotation_degrees,
            )
            try:
                body_img.pack()
            except Exception:
                pass
        except Exception as e:
            print(f"[assemble-rocket] WARNING: failed to apply logo: {e}")

    body_mat = make_image_texture_mat("Mat_Body_Paint", body_img, body_uv_layer, extension="REPEAT")
    assign_material(body, body_mat)

    if write_body_paint_png:
        body_png = Path(out_path).with_suffix("")
        body_png = body_png.parent / f"{body_png.name}_body_paint.png"
        save_image_png(body_img, body_png)
        print(f"[assemble-rocket] Wrote body paint PNG: {body_png}")

fin_min, fin_max = local_aabb(fin_template)
fin_axis_len = fin_max[ai] - fin_min[ai]

# If the fin mesh lies on the chosen spin axis (radial distance ~ 0), rotating won't create
# distinct fins. Print diagnostics and pick a better axis if needed.
rad_x = axis_radial_max(fin_template, "X")
rad_y = axis_radial_max(fin_template, "Y")
rad_z = axis_radial_max(fin_template, "Z")
print(f"[assemble-rocket] fin radial max: X={rad_x:.6f}, Y={rad_y:.6f}, Z={rad_z:.6f}")

spin_axis = fin_spin_axis.upper()
rad_by_axis = {"X": rad_x, "Y": rad_y, "Z": rad_z}
best_axis = max(rad_by_axis.items(), key=lambda kv: kv[1])[0]
if rad_by_axis.get(spin_axis, 0.0) < 1e-5 and best_axis != spin_axis:
    print(
        f"[assemble-rocket] WARNING: fin_spin_axis={spin_axis} has ~zero radial offset; "
        f"rotating around it will stack fins. Falling back to {best_axis}."
    )
    spin_axis = best_axis

# Precompute hinge base point + outward extent from the fin template geometry.
hinge_base = estimate_hinge_point(fin_template, rocket_axis, fin_anchor_t)
radials = [radial_distance_from_axis(v.co, rocket_axis) for v in fin_template.data.vertices]
radial_span = max(radials) - min(radials) if radials else 0.25

# Decide which local axis represents the fin side faces.
if fin_side_axis_mode.upper() == "AUTO":
    fin_side_axis = detect_fin_thickness_axis_by_aabb(fin_template)
else:
    m = fin_side_axis_mode.upper()
    fin_side_axis = {"X": Vector((1.0, 0.0, 0.0)), "Y": Vector((0.0, 1.0, 0.0)), "Z": Vector((0.0, 0.0, 1.0))}[m]

# Create the fin atlas image/material (shared by all fins) if enabled.
fin_atlas_img = None
fin_atlas_mat = None
if use_fin_texture_atlas:
    # 4 cols (fins) x 3 rows (outer/inner/edge)
    cols = 4
    rows = 3
    tile_w = fin_atlas_size // cols
    tile_h = fin_atlas_size // rows

    def _fill_pixels(size, outer_colors, inner_colors, edge_color):
        pixels = [0.0] * (size * size * 4)

        def set_px(x, y, rgba):
            idx = (y * size + x) * 4
            pixels[idx : idx + 4] = list(rgba)

        def fill_tile(col, row, rgba):
            x0 = col * tile_w
            y0 = (rows - 1 - row) * tile_h
            for yy in range(y0, min(y0 + tile_h, size)):
                for xx in range(x0, min(x0 + tile_w, size)):
                    set_px(xx, yy, rgba)

        def draw_fin_index_dots(col, row, dot_count: int):
            """
            Draw 1/2/3 dots in a tile to label Fin_1/2/3.
            We draw a white ring + black center so it's visible on any fin color.
            """
            if dot_count <= 0:
                return
            x0 = col * tile_w
            y0 = (rows - 1 - row) * tile_h
            # Place dots near the top-left of the tile (but inside margin).
            margin_x = max(2, int(tile_w * 0.15))
            margin_y = max(2, int(tile_h * 0.15))
            spacing = max(3, int(tile_w * 0.20))
            r_outer = max(2, int(min(tile_w, tile_h) * 0.075))
            r_inner = max(1, int(r_outer * 0.55))

            for k in range(dot_count):
                cx = x0 + spacing - 0 * dot_count * spacing / 4 + k * spacing
                cy = y0 + tile_h / 2
                draw_marker(pixels, size, cx, cy, r_outer, (0.0, 0.0, 0.0, 1.0))

        for ii in range(4):
            fill_tile(ii, 0, outer_colors[ii])
            fill_tile(ii, 1, inner_colors[ii])
            fill_tile(ii, 2, edge_color)

            # Add visible fin index markers on both sides:
            # - Fin_0: none
            # - Fin_1: 1 dot
            # - Fin_2: 2 dots
            # - Fin_3: 3 dots
            if ii > 0:
                draw_fin_index_dots(ii, 0, ii)  # outer row
                draw_fin_index_dots(ii, 1, ii)  # inner row
        return pixels

    fin_atlas_img = bpy.data.images.get(fin_atlas_image_name)
    if fin_atlas_img is None:
        fin_atlas_img = bpy.data.images.new(
            name=fin_atlas_image_name,
            width=fin_atlas_size,
            height=fin_atlas_size,
            alpha=True,
            float_buffer=True,
        )
    else:
        fin_atlas_img.scale(fin_atlas_size, fin_atlas_size)

    fin_atlas_img.pixels = _fill_pixels(fin_atlas_size, fin_outer_colors, fin_inner_colors, fin_edge_color)
    try:
        fin_atlas_img.pack()
    except Exception:
        pass
    fin_atlas_img.colorspace_settings.name = "sRGB"

    if write_fin_atlas_png:
        atlas_path = Path(out_path).with_suffix("")
        atlas_png = atlas_path.parent / f"{atlas_path.name}_fin_atlas.png"
        atlas_png.parent.mkdir(parents=True, exist_ok=True)
        fin_atlas_img.filepath_raw = str(atlas_png)
        fin_atlas_img.file_format = "PNG"
        try:
            fin_atlas_img.save()
            print(f"[assemble-rocket] Wrote fin atlas PNG: {atlas_png}")
        except Exception as e:
            print(f"[assemble-rocket] WARNING: failed to save fin atlas PNG: {e}")

    fin_atlas_mat = make_image_texture_mat("Mat_Fin_Atlas", fin_atlas_img, fin_uv_layer, extension="CLIP")


def paint_fin_sides_vertex_colors(
    fin_obj,
    local_axis: Vector,
    color_pos_rgba,
    color_neg_rgba,
    edge_rgba,
    layer_name: str = "Col",
    align_threshold: float = 0.92,
):
    """
    "Paint" the fin by writing vertex colors per face:
    - faces whose normal aligns with +local_axis -> color_pos_rgba
    - faces whose normal aligns with -local_axis -> color_neg_rgba
    - all other faces -> edge_rgba
    """
    axis = local_axis.normalized()
    mesh = fin_obj.data

    # Ensure normals are up-to-date.
    mesh.calc_normals()

    vcol = mesh.vertex_colors.get(layer_name)
    if vcol is None:
        vcol = mesh.vertex_colors.new(name=layer_name)
    # Make sure this is the active layer so glTF exporter picks it up (Blender 2.93).
    mesh.vertex_colors.active = vcol
    mesh.vertex_colors.active_index = mesh.vertex_colors.find(vcol.name)
    try:
        mesh.vertex_colors.active_render = vcol
    except Exception:
        pass

    for poly in mesh.polygons:
        d = poly.normal.dot(axis)
        if abs(d) >= align_threshold:
            col = color_pos_rgba if d > 0.0 else color_neg_rgba
        else:
            col = edge_rgba
        for li in poly.loop_indices:
            vcol.data[li].color = col


def fin_tile_uv_rect(fin_index: int, row: int):
    """Return (u0, v0, u1, v1) for a tile in the 4x3 atlas."""
    cols = 4
    rows = 3
    u0 = fin_index / cols
    u1 = (fin_index + 1) / cols
    # v is bottom->top in UV space.
    # Our image fill places row=0 at the TOP, row=rows-1 at the BOTTOM.
    # So the UV rect for row r is:
    # - bottom = (rows-(r+1))/rows
    # - top    = (rows-r)/rows
    v0 = (rows - (row + 1)) / rows
    v1 = (rows - row) / rows
    return (u0, v0, u1, v1)


def assign_fin_uvs_for_atlas(fin_obj, fin_index: int, local_axis: Vector, uv_layer_name: str, align_threshold: float = 0.92):
    """
    Assign UVs so:
    - +axis faces -> outer row (0)
    - -axis faces -> inner row (1)
    - others -> edge row (2)
    """
    axis = local_axis.normalized()
    mesh = fin_obj.data

    if bpy.context.mode != "OBJECT":
        bpy.ops.object.mode_set(mode="OBJECT")

    bm = bmesh.new()
    bm.from_mesh(mesh)
    bm.faces.ensure_lookup_table()
    uv_layer = bm.loops.layers.uv.get(uv_layer_name)
    if uv_layer is None:
        uv_layer = bm.loops.layers.uv.new(uv_layer_name)

    outer_rect = fin_tile_uv_rect(fin_index, 0)
    inner_rect = fin_tile_uv_rect(fin_index, 1)
    edge_rect = fin_tile_uv_rect(fin_index, 2)

    # Use a projection-based UV assignment so it works for tris/quads/ngons.
    # Project onto the plane perpendicular to the side axis and normalize by the fin's bounds.
    side = local_axis.normalized()
    if abs(side.dot(Vector((1.0, 0.0, 0.0)))) > 0.9:
        # side axis ~ X, use YZ plane
        a_idx, b_idx = 1, 2
    elif abs(side.dot(Vector((0.0, 1.0, 0.0)))) > 0.9:
        # side axis ~ Y, use XZ plane
        a_idx, b_idx = 0, 2
    else:
        # side axis ~ Z, use XY plane
        a_idx, b_idx = 0, 1

    verts = [v.co for v in bm.verts]
    a_vals = [v[a_idx] for v in verts] or [0.0]
    b_vals = [v[b_idx] for v in verts] or [0.0]
    a_min, a_max = min(a_vals), max(a_vals)
    b_min, b_max = min(b_vals), max(b_vals)
    a_len = max(1e-6, a_max - a_min)
    b_len = max(1e-6, b_max - b_min)

    def map_into_rect(rect, a_norm: float, b_norm: float):
        u0, v0, u1, v1 = rect
        u = u0 + max(0.0, min(1.0, a_norm)) * (u1 - u0)
        v = v0 + max(0.0, min(1.0, b_norm)) * (v1 - v0)
        return (u, v)

    # Hybrid classification:
    # - Prefer face normal when it's aligned to the side axis (most reliable for the big side faces).
    # - Otherwise fall back to vertex-side "majority vote" (handles base topology where normals get messy).
    s_vals = [v.co.dot(axis) for v in bm.verts] or [0.0]
    s_min, s_max = min(s_vals), max(s_vals)
    s_len = max(1e-6, s_max - s_min)
    # If the mesh straddles 0 along this axis, use 0 as the mid-plane; otherwise use midpoint.
    s_mid = 0.0 if (s_min < 0.0 < s_max) else (s_min + s_max) * 0.5
    s_eps = 0.01 * s_len

    for f in bm.faces:
        d = f.normal.dot(axis)
        if abs(d) >= align_threshold:
            rect = outer_rect if d > 0.0 else inner_rect
        else:
            signed = [loop.vert.co.dot(axis) - s_mid for loop in f.loops]
            pos = sum(1 for s in signed if s > s_eps)
            neg = sum(1 for s in signed if s < -s_eps)
            # Majority vote with a bias toward "side" if one sign dominates.
            if pos > neg * 2 and pos > 0:
                rect = outer_rect
            elif neg > pos * 2 and neg > 0:
                rect = inner_rect
            elif pos > 0 and neg == 0:
                rect = outer_rect
            elif neg > 0 and pos == 0:
                rect = inner_rect
            else:
                rect = edge_rect

        for loop in f.loops:
            co = loop.vert.co
            a_norm = (co[a_idx] - a_min) / a_len
            b_norm = (co[b_idx] - b_min) / b_len
            loop[uv_layer].uv = map_into_rect(rect, a_norm, b_norm)

    bm.to_mesh(mesh)
    bm.free()

# ============================================================
# CREATE FIN PIVOTS + FIN INSTANCES (node instancing)
#
# Goal: keep fins easy to animate via hinge pivots, and keep the body separate.
# Note: restoring per-fin atlas colors requires per-fin UVs, which implies **one fin mesh per fin**
# (4 total). This is still much less duplication than the old "join everything" approach.
# ============================================================

# Root node for fin pivots (matches Bevy name expectations).
bpy.ops.object.empty_add(type="PLAIN_AXES", location=(0.0, 0.0, 0.0))
root_node = bpy.context.object
root_node.name = "Root"
root_node.empty_display_size = 0.05

# Parent body under Root for a clean hierarchy.
body.parent = root_node
body.matrix_parent_inverse = root_node.matrix_world.inverted()

# Optional ghost template: one extra mesh datablock, instanced 4x.
fin_ghost_alpha = 0.12
if use_fin_ghost and use_fin_texture_atlas and fin_atlas_img is not None:
    fin_ghost_template = fin_template.copy()
    fin_ghost_template.data = fin_template.data.copy()
    fin_ghost_template.name = "FinGhostTemplate"
    bpy.context.scene.collection.objects.link(fin_ghost_template)

    ghost_mat = make_image_texture_mat_with_alpha("Mat_Fin_Atlas_Ghost", fin_atlas_img, fin_uv_layer, alpha=fin_ghost_alpha, extension="CLIP")
    assign_material(fin_ghost_template, ghost_mat)
    ensure_uv_layer(fin_ghost_template.data, fin_uv_layer)
    assign_fin_uvs_for_atlas(fin_ghost_template, 0, fin_side_axis, fin_uv_layer, align_threshold=fin_side_align_threshold)

# Build fin pivot empties + instanced fin meshes.
export_objs = [root_node, body]
fins = []
ghost_fins = []

for i in range(fin_count):
    angle = (2 * math.pi / fin_count) * i
    axis_vec = {"X": Vector((1.0, 0.0, 0.0)), "Y": Vector((0.0, 1.0, 0.0)), "Z": Vector((0.0, 0.0, 1.0))}[spin_axis]
    rot_i = Matrix.Rotation(angle, 4, axis_vec)

    # Hinge pivot position and outward axis (same math as old bone placement).
    hinge = rot_i @ hinge_base.copy()
    outward = outward_dir_from_hinge(hinge, rocket_axis)

    # Real fin pivot node: Name matches Bevy's fin lookup ("Root.Fin_i").
    bpy.ops.object.empty_add(type="PLAIN_AXES", location=hinge)
    pivot = bpy.context.object
    pivot.name = f"Root.Fin_{i}"
    pivot.empty_display_size = 0.03
    pivot.parent = root_node
    pivot.matrix_parent_inverse = root_node.matrix_world.inverted()
    # Orient the pivot so animating a single axis behaves like a hinge.
    pivot.rotation_mode = "QUATERNION"
    pivot.rotation_quaternion = hinge_pivot_quat(axis_vec, outward)

    # Fin instance: per-fin mesh data so each fin can have its own atlas UVs / colors.
    fin = fin_template.copy()
    fin.data = fin_template.data.copy()
    fin.name = f"Fin_{i}"
    bpy.context.scene.collection.objects.link(fin)

    if use_fin_texture_atlas and fin_atlas_mat is not None:
        assign_material(fin, fin_atlas_mat)
        ensure_uv_layer(fin.data, fin_uv_layer)
        assign_fin_uvs_for_atlas(fin, i, fin_side_axis, fin_uv_layer, align_threshold=fin_side_align_threshold)

    if place_fins_from_template:
        fin.matrix_world = rot_i @ fin_template.matrix_world
    else:
        # Legacy placement mode (kept for reference / future use).
        fin_axis = body_axis_min + fin_axis_offset_m
        fin.matrix_world = fin_template.matrix_world.copy()
        if rocket_axis.upper() == "X":
            fin.matrix_world.translation.x = fin_axis
        elif rocket_axis.upper() == "Y":
            fin.matrix_world.translation.y = fin_axis
        else:
            fin.matrix_world.translation.z = fin_axis

    # Parent under pivot, preserving world transform (so pivot drives fin motion in-engine).
    bpy.ops.object.select_all(action="DESELECT")
    fin.select_set(True)
    pivot.select_set(True)
    bpy.context.view_layer.objects.active = pivot
    bpy.ops.object.parent_set(type="OBJECT", keep_transform=True)

    export_objs.extend([pivot, fin])
    fins.append(fin)

    # Ghost fin pivot + instance (optional).
    if use_fin_ghost and fin_ghost_template is not None:
        bpy.ops.object.empty_add(type="PLAIN_AXES", location=hinge)
        gpivot = bpy.context.object
        gpivot.name = f"FinGhost_{i}"
        gpivot.empty_display_size = 0.03
        gpivot.parent = root_node
        gpivot.matrix_parent_inverse = root_node.matrix_world.inverted()
        # Same hinge orientation for ghost pivot.
        gpivot.rotation_mode = "QUATERNION"
        gpivot.rotation_quaternion = hinge_pivot_quat(axis_vec, outward)

        g = fin_ghost_template.copy()
        g.data = fin_ghost_template.data  # shared ghost mesh data!
        g.name = f"FinGhostMesh_{i}"
        bpy.context.scene.collection.objects.link(g)
        g.matrix_world = rot_i @ fin_ghost_template.matrix_world
        bpy.ops.object.select_all(action="DESELECT")
        g.select_set(True)
        gpivot.select_set(True)
        bpy.context.view_layer.objects.active = gpivot
        bpy.ops.object.parent_set(type="OBJECT", keep_transform=True)

        export_objs.extend([gpivot, g])
        ghost_fins.append(g)

# ============================================================
# ANCHOR AT BASE: translate all Root children by same amount and apply
# ============================================================
if children_offset_along_axis != 0.0:
    ai = axis_index(rocket_axis)
    for child in list(root_node.children):
        trans = child.location.copy()
        trans[ai] += children_offset_along_axis
        bake_object_translation(child, trans)

# ============================================================
# EXPORT GLB
# ============================================================
bpy.ops.object.select_all(action="DESELECT")
for obj in export_objs:
    try:
        obj.select_set(True)
    except Exception:
        pass
bpy.context.view_layer.objects.active = body

bpy.ops.export_scene.gltf(
    filepath=out_path,
    export_format="GLB",
    use_selection=True,
    export_yup=True,
    export_apply=False,
    export_skins=False,
    export_animations=False,
)

print(f"Exported: {out_path}")
