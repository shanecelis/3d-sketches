use bevy::gltf::GltfAssetLabel;
use bevy::prelude::*;
use bevy_editor_cam::prelude::*;

// const ROCKET_GLB: &str = "blender-rocket.glb";
// const ROCKET_GLB: &str = "assemble-rocket.glb";
const ROCKET_GLB: &str = "assemble-rocket-ghost.glb";
// const ROCKET_GLB: &str = "elodin-rocket.glb";
const FIN_COUNT: usize = 4;

#[derive(Resource, Default)]
struct FinRig {
    fin_entities: [Option<Entity>; FIN_COUNT],
    fin_base_rotations: [Option<Quat>; FIN_COUNT],

    ghost_entities: [Option<Entity>; FIN_COUNT],
    ghost_base_rotations: [Option<Quat>; FIN_COUNT],

    // Per-fin state for lagging the physical fins toward the commanded angle.
    fin_angle: [f32; FIN_COUNT],
}

fn main() {
    App::new()
        .add_plugins(
            DefaultPlugins.set(WindowPlugin {
                primary_window: Some(Window {
                    title: "rocket-wave".to_string(),
                    ..default()
                }),
                ..default()
            }),
        )
        .add_plugins(DefaultEditorCamPlugins)
        .init_resource::<FinRig>()
        .add_systems(Startup, setup)
        .add_systems(Update, (find_fin_nodes, animate_fins))
        .run();
}

fn setup(mut commands: Commands, asset_server: Res<AssetServer>) {
    // Camera
    commands.spawn((
        Camera3d::default(),
        Transform::from_xyz(2.5, 0.0, 1.0).looking_at(Vec3::new(0.0, 0.0, 0.0), Vec3::Y),
        EditorCam::default(),
    ));

    // Light
    commands.spawn((
        DirectionalLight {
            illuminance: 20_000.0,
            shadows_enabled: true,
            ..default()
        },
        Transform::from_xyz(3.0, 8.0, 4.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));

    // Load the first scene from the GLB.
    // Put your file at: assets/rocket.glb
    let rocket_scene: Handle<Scene> = asset_server.load(GltfAssetLabel::Scene(0).from_asset(ROCKET_GLB));
    commands.spawn((Name::new("Rocket"), SceneRoot(rocket_scene)));
}

fn find_fin_nodes(
    mut rig: ResMut<FinRig>,
    named_transforms: Query<(Entity, &Name, &Transform)>,
) {
    for (entity, name, transform) in &named_transforms {
        let n = name.as_str();

        if let Some(i) = fin_index_from_name(n) {
            if rig.fin_entities[i].is_none() {
                rig.fin_entities[i] = Some(entity);
                rig.fin_base_rotations[i] = Some(transform.rotation);
                info_once!("Found fin bone: {} ({entity:?})", n);
            }
        }

        if let Some(i) = ghost_index_from_name(n) {
            if rig.ghost_entities[i].is_none() {
                rig.ghost_entities[i] = Some(entity);
                rig.ghost_base_rotations[i] = Some(transform.rotation);
                info_once!("Found fin ghost bone: {} ({entity:?})", n);
            }
        }
    }
}

fn fin_index_from_name(name: &str) -> Option<usize> {
    // Accept either "Fin_0" or "Root.Fin_0" (and other "<prefix>.Fin_0"-style names).
    let (prefix, idx_str) = name.rsplit_once("Fin_")?;
    if !(prefix.is_empty() || prefix.ends_with('.')) {
        return None;
    }
    let idx: usize = idx_str.parse().ok()?;
    if idx < FIN_COUNT { Some(idx) } else { None }
}

fn ghost_index_from_name(name: &str) -> Option<usize> {
    // Accept exactly "FinGhost_0".."FinGhost_3" (and also allow "<prefix>.FinGhost_0").
    let (prefix, idx_str) = name.rsplit_once("FinGhost_")?;
    if !(prefix.is_empty() || prefix.ends_with('.')) {
        return None;
    }
    let idx: usize = idx_str.parse().ok()?;
    if idx < FIN_COUNT { Some(idx) } else { None }
}

fn animate_fins(time: Res<Time>, mut rig: ResMut<FinRig>, mut transforms: Query<&mut Transform>) {
    // Commanded ("ghost") wave motion. Physical fins lag behind via a simple 1st-order response.
    let t = time.elapsed_secs();
    let dt = time.delta_secs();
    let amplitude = 0.6;
    let speed = 1.6;
    let phase_step = std::f32::consts::TAU / FIN_COUNT as f32;

    for i in 0..FIN_COUNT {
        let commanded = amplitude * (t * speed + i as f32 * phase_step).sin();

        // If ghost bones exist, animate them directly as the commanded target.
        if let (Some(e), Some(base)) = (rig.ghost_entities[i], rig.ghost_base_rotations[i]) {
            if let Ok(mut tr) = transforms.get_mut(e) {
                tr.rotation = base * Quat::from_rotation_y(commanded);
            }
        }

        // Physical fin: if ghost exists, lag toward commanded; otherwise just use commanded.
        let target = if rig.ghost_entities[i].is_some() {
            commanded
        } else {
            commanded
        };

        // First-order lag: dθ/dt = (target-θ)/tau
        // tau controls how "sluggish" the fin is.
        let tau = 0.48_f32;
        let alpha = if tau > 0.0 { (dt / tau).min(1.0) } else { 1.0 };
        rig.fin_angle[i] = rig.fin_angle[i] + (target - rig.fin_angle[i]) * alpha;

        if let (Some(e), Some(base)) = (rig.fin_entities[i], rig.fin_base_rotations[i]) {
            if let Ok(mut tr) = transforms.get_mut(e) {
                tr.rotation = base * Quat::from_rotation_y(rig.fin_angle[i]);
            }
        }
    }
}
