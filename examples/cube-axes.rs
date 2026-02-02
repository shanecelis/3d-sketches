use bevy::gltf::GltfAssetLabel;
use bevy::prelude::*;
use bevy::picking::prelude::{Click, Pointer};

fn main() {
    App::new()
        // MeshPickingPlugin is not a default plugin.
        .add_plugins((DefaultPlugins, MeshPickingPlugin))
        .add_plugins(bevy_inspector_egui::bevy_egui::EguiPlugin::default())
        .add_plugins(bevy_inspector_egui::quick::WorldInspectorPlugin::new())
        .add_systems(Startup, setup)
        // glTF scenes spawn meshes async; attach click observers as they appear.
        .add_systems(Update, attach_click_observers_to_spawned_meshes)
        .run();
}

#[derive(Component)]
struct AxesCubeRoot;

fn setup(mut commands: Commands, asset_server: Res<AssetServer>) {
    commands.spawn((
        Camera3d::default(),
        Transform::from_xyz(1.8, 1.4, 2.6).looking_at(Vec3::ZERO, Vec3::Y),
    ));

    commands.spawn((
        Camera2d::default(),
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
    let cube = commands
        .spawn((
            Name::new("AxesCube"),
            AxesCubeRoot,
            // SceneRoot(scene),
            InheritedVisibility::default(),
            Transform::default(),
        ))
        .id();

    // Face labels: 6 Text2d entities as children so they rotate with the cube.
    //
    // Coordinate convention:
    // - Front: +Y
    // - Right: +X
    // - Top: +Z
    // - Back: -Y
    // - Left: -X
    // - Bottom: -Z
    let d = 0.62; // slightly in front of the cube face (cube_size=1 => half=0.5)
    // let s = 0.20; // text scale
    let s = 1.00; // text scale
    // Basic usage.
    commands.spawn(Text2d::new("hello world!"));

    commands.entity(cube).with_children(|p| {
        p.spawn((
            Name::new("Label_Front"),
            Text2d::new("front"),
            TextFont::default(),
            Transform::from_translation(Vec3::new(0.0, d, 0.0))
                .with_rotation(Quat::from_rotation_x(-std::f32::consts::FRAC_PI_2))
                .with_scale(Vec3::splat(s)),
        ));
        p.spawn((
            Name::new("Label_Right"),
            Text2d::new("right"),
            Transform::from_translation(Vec3::new(d, 0.0, 0.0))
                .with_rotation(Quat::from_rotation_y(std::f32::consts::FRAC_PI_2))
                .with_scale(Vec3::splat(s)),
        ));
        p.spawn((
            Name::new("Label_Top"),
            Text2d::new("top"),
            Transform::from_translation(Vec3::new(0.0, 0.0, d)).with_scale(Vec3::splat(s)),
        ));
        p.spawn((
            Name::new("Label_Back"),
            Text2d::new("back"),
            Transform::from_translation(Vec3::new(0.0, -d, 0.0))
                .with_rotation(
                    Quat::from_rotation_x(std::f32::consts::FRAC_PI_2)
                        * Quat::from_rotation_z(std::f32::consts::PI),
                )
                .with_scale(Vec3::splat(s)),
        ));
        p.spawn((
            Name::new("Label_Left"),
            Text2d::new("left"),
            Transform::from_translation(Vec3::new(-d, 0.0, 0.0))
                .with_rotation(Quat::from_rotation_y(-std::f32::consts::FRAC_PI_2))
                .with_scale(Vec3::splat(s)),
        ));
        p.spawn((
            Name::new("Label_Bottom"),
            Text2d::new("bottom"),
            Transform::from_translation(Vec3::new(0.0, 0.0, -d))
                .with_rotation(Quat::from_rotation_y(std::f32::consts::PI))
                .with_scale(Vec3::splat(s)),
        ));
    });
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

