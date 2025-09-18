# Synthetic MPR

> For users: TM.txt is incomplete, it's still runable but the results will be non-physical. Please contact Xutao Xu for a realistic TM.txt if you want.

## A quick start

1. Specify `Target`, `Aperture`, `Magnets` and `Focalplane`, they are all **class** objects defined in `components.py`.
2. Run in the following order: `Beam.generate()` (or other alternatives), `Beam.trans()` and `Beam.hit()`. All methods above are mentioned under **Beam: class** in `processes.py`.
3. Do the following visualization and statistics as you wish!

```python
# There are two ways to specify the geometry of one target: 
# 1) by giving geometry directly; 2) by giving area and H_W_ratio.
target1 = Target(type = 'CH2', mass_thickness = 1.0, shape = 'rectangle', geometry = [width, height])
target2 = Target(type = 'CH2', mass_thickness = 1.0, shape = 'rectangle', area = 1.0, H_W_ratio = 3/2)

# Same situation when it comes to the specification of an aperture:
aperture = Aperture(distance = 1.0, shape = 'circle', area = 1.0, H_W_ratio = 1.0)

magnets = Magnets(file_path = "TM.txt", reference_energy = 10)

# Two kinds of focal plane: normal (perpendicular to the reference beam, straight); arbitrary.
focalplane1 = FocalPlane(type = 'normal', position = 1.0)
focalplane2 = FocalPlane(type = 'arbitrary', geomotry = geometry)

# Beam generation - using Target and Aperture
beam1 = Beam.generate(Target, Aperture, energy = 10, Npart = 100, cross_section = cross_section)
beam2 = Beam.generate_monoenergetic_parrallel(Target, energy = 10, Npart = 100)
beam3 = Beam.generate_monoenergetic_cone(energy = 10, Npart = 100, a_max = 0.01)

beam_tranported = beam1.trans(magnets)

record = beam_transported.hit(focalplane1)
```

## Schematic Diagram

For you to better understand the setting of this synthetic system, see the schematic diagram down below.

![4 components and 3 processes, that's all you could expect!](./readme_imgs/Diagram.png)

## Useful integrated tools in `mpr_system.py`

Simply using `Target`, `Aperture`, `Magnets` and `Focalplane`, you can specify a new class object: `MPR`. It's a great way to represent your MPR system!

```python
mpr = MPR(target, aperture, magnet, focalplane)
```

Or, if you're not sure about where to place your focal plane, it's OK to initialize without a `FocalPlane` object:

```python
mpr = MPR(target, aperture, magnet)
```

Up to now, the `MPR` object has the following functions:

1. `MPR.performance()`: calculates l-E relation and energy resolution in a given energy range.
2. `MPR.optimal_focalplane()`: finds optimal **straight** focal plane position (including position and tilt angle) under given boundaries. This method would automatically **overwrite** self.focalplane. So, if you initialized `MPR` without providing `FocalPlane`, you should run this method before others.
3. `MPR.save_focalplane()`: saves the current focal plane's geometry (list of (x, z) tuples) to the given path.
4. `MPR.response_matrix()`: calculates response matrix of the MPR system.

```python
mpr.performance()

# The merit, which is minimized during optimization, can be difined in two ways:
# 1) uniformly average of energy resolution at different energy; 2) weighted averge of ~  
mpr.optimal_focalplane(merit_weight = 'uniform', fp_position = np.linspace(0, 0.5, 11), fp_angle = np.linspace(10, 30, 6))
mpr.optimal_focalplane(merit_weight = 'whatever', fp_position = np.linspace(0, 0.5, 11), fp_angle = np.linspace(10, 30, 6), merit_weight_lst = [1,2,3,4,5,6,5,4,3,2,1])

mpr.save_focalplane(save_path = save_path)

mpr.response_matrix(save_path = save_path, plot_save = True)
```
