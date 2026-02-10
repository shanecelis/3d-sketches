BLENDER ?= /Applications/Blender.app/Contents/MacOS/Blender

all: assets/blender-rocket.glb assets/axes-cube.glb assets/assemble-rocket.glb assets/crazyflie-rig.glb assets/FinGhost.glb assets/assemble-rocket-ghost.glb

# Build any `assets/<name>.glb` by running `<name>.py` and passing the output path after `--`.
assets/%.glb: %.py | assets
	$(BLENDER) --background --factory-startup --python $< -- $@

assets/assemble-rocket-ghost.glb: assemble-rocket.py
	$(BLENDER) --background --factory-startup --python $< -- $@ Body.glb Fin.glb --add-ghost-fins

assets:
	mkdir -p $@
