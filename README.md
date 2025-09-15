# Synthetic MPR

> For users: TM.txt is incomplete, it's still runable but the results will be non-physical. Please contact Xutao Xu for a realistic TM.txt if you want.

## A quick start

1. Specify `Target`, `Aperture`, `Magnets` and `Focalplane`, they are all **class** objects defined in `components.py`.
2. Run in the following order: `Beam.generate()` (or other alternatives), `Beam.trans()` and `Beam.hit()`. All methods above are mentioned under **Beam: class** in `processes.py`.
3. Do the following visualization and statistics as you wish!

## Schematic Diagram

For you to better understand the setting of this synthetic system, see the schematic diagram down below.

![4 components and 3 processes, that's all you could expect!](./readme_imgs/Diagram.png)

## Useful integrated tools in `mpr_system.py`

Simply using `Target`, `Aperture`, `Magnets` and `Focalplane`, you can specify a new class object: `MPR`. It's a great way to represent your MPR system!

Up to now, it has the following functions:

1. `MPR.performance()`: calculates l-E relation and energy resolution in a given energy range.
2. `MPR.optimal_focalplane()`: finds optimal **straight** focal plane position (including position and tilt angle) under given boundaries.
