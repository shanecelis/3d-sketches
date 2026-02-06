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
#   blender ... --python assemble-rocket.py -- out.glb Body.glb Fin.glb
argv = sys.argv
if "--" in argv:
    user_args = argv[argv.index("--") + 1 :]
else:
    user_args = []

script_dir = Path(__file__).resolve().parent
body_glb_path = Path(user_args[1]).expanduser().resolve() if len(user_args) >= 2 else (script_dir / "Body.glb")
fin_glb_path = Path(user_args[2]).expanduser().resolve() if len(user_args) >= 3 else (script_dir / "Fin.glb")

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


def make_image_texture_mat(name: str, image: bpy.types.Image, uv_layer_name: str):
    """Material that displays an Image Texture using the given UV map (reliable Blender preview)."""
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
    tex.extension = "CLIP"

    uv = nt.nodes.new("ShaderNodeUVMap")
    uv.location = (-260, -40)
    uv.uv_map = uv_layer_name

    nt.links.new(uv.outputs["UV"], tex.inputs["Vector"])
    nt.links.new(tex.outputs["Color"], bsdf.inputs["Base Color"])
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])
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


def assign_cylindrical_uvs(mesh_obj, axis: str, uv_layer_name: str, u_offset_degrees: float = 0.0):
    """
    Assign simple cylindrical UVs:
    - U wraps around the rocket axis
    - V runs along the rocket axis from min..max
    """
    axis = axis.upper()
    mesh = mesh_obj.data
    if bpy.context.mode != "OBJECT":
        bpy.ops.object.mode_set(mode="OBJECT")

    bm = bmesh.new()
    bm.from_mesh(mesh)
    bm.faces.ensure_lookup_table()
    uv_layer = bm.loops.layers.uv.get(uv_layer_name)
    if uv_layer is None:
        uv_layer = bm.loops.layers.uv.new(uv_layer_name)

    vmin, vmax = local_aabb(mesh_obj)
    if axis == "X":
        amin, amax = vmin.x, vmax.x
    elif axis == "Y":
        amin, amax = vmin.y, vmax.y
    else:
        amin, amax = vmin.z, vmax.z
    alen = max(1e-6, amax - amin)

    two_pi = 2.0 * math.pi
    u_offset = (u_offset_degrees / 360.0) % 1.0
    for f in bm.faces:
        for loop in f.loops:
            co = loop.vert.co
            if axis == "X":
                u = (math.atan2(co.z, co.y) / two_pi) + 0.5
                v = (co.x - amin) / alen
            elif axis == "Y":
                u = (math.atan2(co.x, co.z) / two_pi) + 0.5
                v = (co.y - amin) / alen
            else:  # "Z"
                u = (math.atan2(co.y, co.x) / two_pi) + 0.5
                v = (co.z - amin) / alen
            u = (u + u_offset) % 1.0
            loop[uv_layer].uv = (u, v)

    bm.to_mesh(mesh)
    bm.free()


def create_body_paint_image(name: str, size: int, band_v0: float, band_v1: float, sector_colors):
    """
    Create an image where:
    - everything is black
    - band region (v in [band_v0..band_v1]) is split into 4 U sectors with colors
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

    black = (0.0, 0.0, 0.0, 1.0)

    # Fill all black first.
    for y in range(size):
        for x in range(size):
            set_px(x, y, black)

    y0 = max(0, min(size, int(round(size * max(0.0, min(1.0, band_v0))))))
    y1 = max(0, min(size, int(round(size * max(0.0, min(1.0, band_v1))))))
    if y1 < y0:
        y0, y1 = y1, y0
    if y1 == y0:
        y1 = min(size, y0 + 1)

    for y in range(y0, y1):
        for x in range(size):
            sector = min(3, int((x / size) * 4))
            set_px(x, y, sector_colors[sector])

    img.pixels = pixels
    try:
        img.pack()
    except Exception:
        pass
    img.colorspace_settings.name = "sRGB"
    return img


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

# Body paint via texture (black body + colored band near the tail).
use_body_paint_texture = True
body_paint_image_name = "BodyPaint"
body_paint_size = 1024  # px (square)
# Band is measured from the body's "bottom" along rocket_axis (min axis value).
body_band_start_m = 0.0
body_band_length_m = 0.0508  # 2 inches
# Rotate the colored sectors around the rocket axis (degrees).
# Positive values rotate the pattern in the same direction as a positive rotation about rocket_axis.
body_sector_u_offset_degrees = 135.0
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

# Rotate the body to the desired orientation and bake it in.
bake_object_rotation(body, (0.0, math.radians(body_y_rotation_degrees), 0.0))

# If you want to preserve imported body materials, comment out the next line.
assign_material(body, mat_body)

# Ensure the active object's UV map exists before we join meshes later.
# (Join tends to preserve UV maps from the active object; we keep names consistent.)
ensure_uv_layer(body.data, fin_uv_layer)

body_min, body_max = local_aabb(body)
ai = axis_index(rocket_axis)
body_axis_min = body_min[ai]
body_axis_max = body_max[ai]

# Body paint: black base + colored band near tail (matches fin colors).
if use_body_paint_texture:
    assign_cylindrical_uvs(body, rocket_axis, fin_uv_layer, u_offset_degrees=body_sector_u_offset_degrees)

    body_axis_len = max(1e-6, body_axis_max - body_axis_min)
    band_v0 = (body_band_start_m) / body_axis_len
    band_v1 = (body_band_start_m + body_band_length_m) / body_axis_len
    # Use the fin outer colors as the 4 sectors around the rocket.
    body_img = create_body_paint_image(body_paint_image_name, body_paint_size, band_v0, band_v1, fin_outer_colors)
    body_mat = make_image_texture_mat("Mat_Body_Paint", body_img, fin_uv_layer)
    assign_material(body, body_mat)

    if write_body_paint_png:
        body_png = Path(out_path).with_suffix("")
        body_png = body_png.parent / f"{body_png.name}_body_paint.png"
        body_png.parent.mkdir(parents=True, exist_ok=True)
        body_img.filepath_raw = str(body_png)
        body_img.file_format = "PNG"
        try:
            body_img.save()
            print(f"[assemble-rocket] Wrote body paint PNG: {body_png}")
        except Exception as e:
            print(f"[assemble-rocket] WARNING: failed to save body paint PNG: {e}")

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

        for ii in range(4):
            fill_tile(ii, 0, outer_colors[ii])
            fill_tile(ii, 1, inner_colors[ii])
            fill_tile(ii, 2, edge_color)
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

    fin_atlas_mat = make_image_texture_mat("Mat_Fin_Atlas", fin_atlas_img, fin_uv_layer)


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

    def rect_uvs(rect):
        u0, v0, u1, v1 = rect
        return [Vector((u0, v0)), Vector((u1, v0)), Vector((u1, v1)), Vector((u0, v1))]

    for f in bm.faces:
        d = f.normal.dot(axis)
        if abs(d) >= align_threshold:
            rect = outer_rect if d > 0.0 else inner_rect
        else:
            rect = edge_rect

        uvs = rect_uvs(rect)
        loops = list(f.loops)
        for li, loop in enumerate(loops):
            loop[uv_layer].uv = uvs[li % 4]

    bm.to_mesh(mesh)
    bm.free()

# ============================================================
# CREATE FINS (distributed around rocket)
# Each fin is fully weighted to its corresponding fin bone vertex group.
#
# Assumptions about `Fin.glb` orientation (recommended):
# - local +X: outward direction
# - local +Y: thickness
# - local +Z: rocket axis / fin height
# - origin: near the hinge/pivot (optional but helpful)
# ============================================================
fins = []
for i in range(fin_count):
    angle = (2 * math.pi / fin_count) * i
    fin = fin_template.copy()
    fin.data = fin_template.data.copy()
    fin.name = f"FinObj_{i}"
    bpy.context.scene.collection.objects.link(fin)

    # Color the fin sides.
    # Convention: fin "sides" are the two broad faces whose normals align with ±fin_side_axis.
    if use_fin_texture_atlas:
        assign_material(fin, fin_atlas_mat)
        ensure_uv_layer(fin.data, fin_uv_layer)
        assign_fin_uvs_for_atlas(fin, i, fin_side_axis, fin_uv_layer, align_threshold=fin_side_align_threshold)
    elif use_vertex_colors_for_fin_sides:
        # Use a single material so vertex colors show "as-is" in most engines.
        assign_material(fin, make_vertex_color_mat("Mat_Fin_VCol", fin_vertex_color_layer))
        paint_fin_sides_vertex_colors(
            fin,
            fin_side_axis,
            fin_outer_colors[i],
            fin_inner_colors[i],
            fin_edge_color,
            layer_name=fin_vertex_color_layer,
        )
    else:
        # Material-based split (older approach).
        mat_fin_outer = make_mat(f"Mat_Fin{i}_Outer", fin_outer_colors[i])
        mat_fin_inner = make_mat(f"Mat_Fin{i}_Inner", fin_inner_colors[i])
        assign_fin_two_sided_materials(
            fin,
            Vector((0.0, 1.0, 0.0)),
            mat_fin_outer,
            mat_fin_inner,
            mat_fin_edge,
        )

    if place_fins_from_template:
        # Fin.glb is already placed correctly; just duplicate and spin copies around the rocket axis
        # *in world space* so any imported offset from the axis is preserved.
        fin.matrix_world = Matrix.Rotation(angle, 4, spin_axis) @ fin_template.matrix_world

        bpy.ops.object.select_all(action="DESELECT")
        fin.select_set(True)
        bpy.context.view_layer.objects.active = fin
        bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)
    else:
        # Legacy placement mode (kept for reference / future use).
        fin_axis = body_axis_min + fin_axis_offset_m
        if rocket_axis.upper() == "X":
            fin.location.x = fin_axis
        elif rocket_axis.upper() == "Y":
            fin.location.y = fin_axis
        else:
            fin.location.z = fin_axis

        bpy.ops.object.select_all(action="DESELECT")
        fin.select_set(True)
        bpy.context.view_layer.objects.active = fin
        bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)

    # Match the bone name we create below (so skinning works in glTF).
    weight_all_verts(fin, f"Root.Fin_{i}", 1.0)
    fins.append(fin)

# ============================================================
# JOIN INTO ONE MESH
# ============================================================
for obj in [body] + fins:
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
if rocket_axis.upper() == "X":
    root_bone.head = (body_axis_min, 0, 0)
    root_bone.tail = (body_axis_max, 0, 0)
elif rocket_axis.upper() == "Y":
    root_bone.head = (0, body_axis_min, 0)
    root_bone.tail = (0, body_axis_max, 0)
else:
    root_bone.head = (0, 0, body_axis_min)
    root_bone.tail = (0, 0, body_axis_max)

# Fin bones: hinge at (near) fin base, tail outward
for i in range(fin_count):
    angle = (2 * math.pi / fin_count) * i

    # Fin.glb provides the base fin placement. We estimate the hinge from fin geometry,
    # then rotate that hinge point around the rocket axis for each fin.
    hinge = hinge_base.copy()
    hinge = Matrix.Rotation(angle, 4, spin_axis) @ hinge

    outward = outward_dir_from_hinge(hinge, rocket_axis)

    # Name bones to match the requested convention in-engine.
    b = arm.edit_bones.new(f"Root.Fin_{i}")
    b.head = hinge
    b.tail = hinge + outward * (radial_span * 0.6)
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
# BAKE A GLOBAL ROTATION (so engines don't need to rotate it)
# Rotate the whole model +90° around Y and APPLY the transform.
# ============================================================
def bake_global_y_rotation(mesh_obj, armature_obj, angle_degrees: float):
    """
    Bake a world-space rotation into both the mesh object data and the armature,
    while preserving skinning/parenting.
    """
    angle = math.radians(angle_degrees)
    rot = Matrix.Rotation(angle, 4, "Y")

    # Ensure we're in OBJECT mode.
    if bpy.context.mode != "OBJECT":
        bpy.ops.object.mode_set(mode="OBJECT")

    # Temporarily unparent the mesh so applying transforms on the armature doesn't
    # implicitly compensate the child and cancel out the intended bake.
    bpy.ops.object.select_all(action="DESELECT")
    mesh_obj.select_set(True)
    bpy.context.view_layer.objects.active = mesh_obj
    if mesh_obj.parent is not None:
        bpy.ops.object.parent_clear(type="CLEAR_KEEP_TRANSFORM")

    # Rotate both objects in world space.
    mesh_obj.matrix_world = rot @ mesh_obj.matrix_world
    armature_obj.matrix_world = rot @ armature_obj.matrix_world

    # Apply rotation to bake into data.
    for obj in (mesh_obj, armature_obj):
        bpy.ops.object.select_all(action="DESELECT")
        obj.select_set(True)
        bpy.context.view_layer.objects.active = obj
        bpy.ops.object.transform_apply(location=False, rotation=True, scale=False)

    # Re-parent mesh back to armature, preserving world transform.
    mesh_obj.parent = armature_obj
    mesh_obj.matrix_parent_inverse = armature_obj.matrix_world.inverted()


# Rotate the exported rocket so it faces the desired direction in-engine.
bake_global_y_rotation(rocket, rig, 180.0)

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
    export_colors=True,
    export_skins=True,
    export_animations=False,
)

print(f"Exported: {out_path}")
