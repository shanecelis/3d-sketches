# 3d-sketches

Small Bevy + Blender playground for procedural models, rigging, and picking.

It's more of a sketchpad than a proper project for external use.

## What’s in here

- **Bevy app** (`src/main.rs`)
  - Loads `assets/rocket.glb`
  - Finds fin bones by name (`Fin_0..Fin_3` or `Root.Fin_0..Root.Fin_3`)
  - Animates those bones procedurally.

- **Bevy examples** (`examples/`)
  - `cube-axes.rs`: loads the axes widget and uses Bevy picking to print clicked mesh `Name`s
  - `crazyflie-spin.rs`: loads the rigged Crazyflie and spins prop bones

- **Blender scripts**
  - `blender-rocket.py`: generates a simple rocket with an armature and fin bones, exports GLB
  - `assemble-rocket.py`: assembles a rocket from `Body.glb` + `Fin.glb`, rigs fins, and bakes textures
    - Generates inspection textures in `assets/` (body paint, fin atlas)
  - `axes-cube.py`: generates an “axes cube” widget GLB with separately pickable faces/borders/corners
  - `crazyflie-rig.py`: imports `crazyflie.glb`, adds an armature for `propeller*` meshes, exports `crazyflie-rigged.glb`

## Assets

`assets/` contains exported `.glb` files and generated textures used for inspection.

- `assets/rocket.glb`: runtime rocket model for the main Bevy app
- `assets/axes-cube.glb`: axes widget
- `assets/assemble-rocket.glb`: assembled/painted rocket output
- `assets/crazyflie-rigged.glb`: Crazyflie with prop armature

## Common commands

- **Run main Bevy app**

```bash
cargo run
```

- **Run an example**

```bash
cargo run --example cube-axes
cargo run --example crazyflie-spin
```

- **Generate GLBs via Blender** (macOS default path; see `Makefile`)

```bash
make assets/axes-cube.glb
make assets/blender-rocket.glb
make assets/assemble-rocket.glb
make assets/crazyflie-rigged.glb
```

