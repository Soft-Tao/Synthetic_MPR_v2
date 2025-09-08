# Synthetic MPR

> For users: TM.txt is incomplete. Please contact Xutao Xu.

## A quick start

1. Specify `Target`, `Aperture`, `Magnets` and `Focalplane`, they are all **class** objects defined in `components.py`.
2. Run in the following order: `Beam.generate()` (or other alternatives), `Beam.trans()` and `Beam.hit()`. All methods above are mentioned under **Beam: class** in `processes.py`.
3. Do the following visualization and statistics as you wish!

## Schematic Diagram

For you to better understand the setting of this synthetic system, see the schematic diagram down below.

![4 components and 3 processes, that's all you could expect!](./readme_imgs/Diagram.png)
