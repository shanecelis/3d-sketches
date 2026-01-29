use bevy::gltf::GltfAssetLabel;
use bevy::prelude::*;
use bevy_editor_cam::prelude::*;

const ROCKET_GLB: &str = "rocket.glb";
// const ROCKET_GLB: &str = "elodin-rocket.glb";
const FIN_COUNT: usize = 4;

#[derive(Resource, Default)]
struct FinBones {
    entities: [Option<Entity>; FIN_COUNT],
    base_rotations: [Option<Quat>; FIN_COUNT],
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
        .init_resource::<FinBones>()
        .add_systems(Startup, setup)
        .add_systems(Update, (find_fin_bones, animate_fins))
        .run();
}

fn setup(mut commands: Commands, asset_server: Res<AssetServer>) {
    // Camera
    commands.spawn((
        Camera3d::default(),
        Transform::from_xyz(0.0, 1.5, 6.0).looking_at(Vec3::new(0.0, 1.0, 0.0), Vec3::Y),
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

fn find_fin_bones(
    mut fins: ResMut<FinBones>,
    named_transforms: Query<(Entity, &Name, &Transform)>,
) {
    if fins.entities.iter().all(|e| e.is_some()) {
        return;
    }

    for (entity, name, transform) in &named_transforms {
        let Some(i) = fin_index_from_name(name.as_str()) else {
            continue;
        };
        if fins.entities[i].is_some() {
            continue;
        }
        fins.entities[i] = Some(entity);
        fins.base_rotations[i] = Some(transform.rotation);
        info_once!("Found fin bone: {} ({entity:?})", name.as_str());
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

fn animate_fins(time: Res<Time>, fins: Res<FinBones>, mut transforms: Query<&mut Transform>) {
    // Simple "wave" motion. If the axis is wrong for your model, swap Z for X/Y.
    let t = time.elapsed_secs();
    let amplitude = 0.6;
    let speed = 1.6;
    let phase_step = std::f32::consts::TAU / FIN_COUNT as f32;

    for i in 0..FIN_COUNT {
        let (Some(entity), Some(base_rotation)) = (fins.entities[i], fins.base_rotations[i]) else {
            continue;
        };
        let Ok(mut fin_transform) = transforms.get_mut(entity) else {
            continue;
        };

        let angle = amplitude * (t * speed + i as f32 * phase_step).sin();
        fin_transform.rotation = base_rotation * Quat::from_rotation_y(angle);
    }
}
