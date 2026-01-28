use bevy::gltf::GltfAssetLabel;
use bevy::prelude::*;

const ROCKET_GLB: &str = "blender-rocket.glb";
const FIN_BONE_NAME: &str = "Root.Fin_0";

#[derive(Resource, Default)]
struct FinBone {
    entity: Option<Entity>,
    base_rotation: Option<Quat>,
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
        .init_resource::<FinBone>()
        .add_systems(Startup, setup)
        .add_systems(Update, (find_fin_bone, animate_fin))
        .run();
}

fn setup(mut commands: Commands, asset_server: Res<AssetServer>) {
    // Camera
    commands.spawn((
        Camera3d::default(),
        Transform::from_xyz(0.0, 1.5, 6.0).looking_at(Vec3::new(0.0, 1.0, 0.0), Vec3::Y),
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

fn find_fin_bone(
    mut fin: ResMut<FinBone>,
    named_transforms: Query<(Entity, &Name, &Transform)>,
) {
    if fin.entity.is_some() {
        return;
    }

    for (entity, name, transform) in &named_transforms {
        if name.as_str() == FIN_BONE_NAME {
            fin.entity = Some(entity);
            fin.base_rotation = Some(transform.rotation);
            info_once!("Found fin bone: {FIN_BONE_NAME} ({entity:?})");
            break;
        }
    }
}

fn animate_fin(time: Res<Time>, fin: Res<FinBone>, mut transforms: Query<&mut Transform>) {
    let (Some(entity), Some(base_rotation)) = (fin.entity, fin.base_rotation) else {
        return;
    };

    let Ok(mut t) = transforms.get_mut(entity) else {
        return;
    };

    // Simple "wave" motion. If the axis is wrong for your model, swap Z for X/Y.
    let angle = 0.6 * time.elapsed_secs().sin();
    t.rotation = base_rotation * Quat::from_rotation_z(angle);
}
