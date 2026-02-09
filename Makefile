BLENDER ?= /Applications/Blender.app/Contents/MacOS/Blender

all: assets/blender-rocket.glb assets/axes-cube.glb assets/assemble-rocket.glb assets/crazyflie-rig.glb assets/crazyflie-rigged.glb

# Rig crazyflie props and output a stable name.
assets/crazyflie-rigged.glb: crazyflie-rig.py crazyflie.glb | assets
	$(BLENDER) --background --factory-startup --python crazyflie-rig.py -- $@ crazyflie.glb

# Build any `assets/<name>.glb` by running `<name>.py` and passing the output path after `--`.
assets/%.glb: %.py | assets
	$(BLENDER) --background --factory-startup --python $< -- $@

assets:
	mkdir -p $@
