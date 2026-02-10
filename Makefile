BLENDER ?= /Applications/Blender.app/Contents/MacOS/Blender

all: assets/blender-rocket.glb assets/axes-cube.glb assets/assemble-rocket.glb assets/crazyflie-rig.glb


# Build any `assets/<name>.glb` by running `<name>.py` and passing the output path after `--`.
assets/%.glb: %.py | assets
	$(BLENDER) --background --factory-startup --python $< -- $@

assets:
	mkdir -p $@
