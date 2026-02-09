use bevy::gltf::GltfAssetLabel;
use bevy::prelude::*;

const SCENE: &str = "crazyflie-rigged.glb";

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_systems(Startup, setup)
        .add_systems(Update, (find_props, spin_props))
        .run();
}

#[derive(Resource, Default)]
struct Props {
    entities: Vec<Entity>,
    base_rot: Vec<Quat>,
}

fn setup(mut commands: Commands, asset_server: Res<AssetServer>) {
    commands.spawn((
        Camera3d::default(),
        Transform::from_xyz(0.3, 0.2, 0.6).looking_at(Vec3::ZERO, Vec3::Y),
    ));

    commands.spawn((
        DirectionalLight {
            illuminance: 30_000.0,
            shadows_enabled: true,
            ..default()
        },
        Transform::from_xyz(1.0, 2.0, 1.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));

    let scene: Handle<Scene> = asset_server.load(GltfAssetLabel::Scene(0).from_asset(SCENE));
    commands.spawn((Name::new("Crazyflie"), SceneRoot(scene)));

    commands.init_resource::<Props>();
}

fn find_props(mut props: ResMut<Props>, q: Query<(Entity, &Name, &Transform), Added<Name>>) {
    // Pick up any spawned propeller nodes/bones.
    for (e, name, t) in &q {
        let n = name.as_str().to_lowercase();
        if n.starts_with("root.propeller_") || n.starts_with("propeller") {
            props.entities.push(e);
            props.base_rot.push(t.rotation);
            info!("Found prop node: {} ({e:?})", name.as_str());
        }
    }
}

fn spin_props(time: Res<Time>, props: Res<Props>, mut transforms: Query<&mut Transform>) {
    let w = 40.0; // rad/s
    for (i, &e) in props.entities.iter().enumerate() {
        let Ok(mut t) = transforms.get_mut(e) else { continue };
        let base = props.base_rot.get(i).copied().unwrap_or(Quat::IDENTITY);
        // Crazyflie rig script uses local +Y axis for spin.
        t.rotation = base * Quat::from_rotation_y(time.elapsed_secs() * w);
    }
}

