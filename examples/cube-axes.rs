use bevy::gltf::GltfAssetLabel;
use bevy::prelude::*;
use bevy::picking::prelude::{Click, Pointer};

fn main() {
    App::new()
        // MeshPickingPlugin is not a default plugin.
        .add_plugins((DefaultPlugins, MeshPickingPlugin))
        .add_systems(Startup, setup)
        // glTF scenes spawn meshes async; attach click observers as they appear.
        .add_systems(Update, attach_click_observers_to_spawned_meshes)
        .run();
}

fn setup(mut commands: Commands, asset_server: Res<AssetServer>) {
    commands.spawn((
        Camera3d::default(),
        Transform::from_xyz(1.8, 1.4, 2.6).looking_at(Vec3::ZERO, Vec3::Y),
    ));

    commands.spawn((
        DirectionalLight {
            illuminance: 25_000.0,
            shadows_enabled: true,
            ..default()
        },
        Transform::from_xyz(2.0, 4.0, 2.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));

    let scene: Handle<Scene> = asset_server.load(GltfAssetLabel::Scene(0).from_asset("axes-cube.glb"));
    commands.spawn((Name::new("AxesCube"), SceneRoot(scene)));
}

#[derive(Component)]
struct ClickObserved;

/// glTF scenes spawn many entities; attach a click observer to each spawned mesh entity.
///
/// With `MeshPickingPlugin`, entities with `Mesh3d` are pickable by default.
fn attach_click_observers_to_spawned_meshes(
    mut commands: Commands,
    meshes: Query<Entity, (Added<Mesh3d>, Without<ClickObserved>)>,
) {
    for e in &meshes {
        commands
            .entity(e)
            .insert(ClickObserved)
            .observe(print_clicked_name);
    }
}

fn print_clicked_name(click: On<Pointer<Click>>, names: Query<&Name>) {
    if let Ok(name) = names.get(click.entity) {
        info!("Clicked: {}", name.as_str());
    } else {
        info!("Clicked entity {:?}", click.entity);
    }
}

