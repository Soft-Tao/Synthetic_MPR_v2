Overview
This project simulates particle beam behavior through a target, aperture, and magnetic transport system, with visualization of particle impacts on a focal plane. The code consists of three main components:

components.py: Physical component definitions
processes.py: Core simulation processes
syn_MPR.py: Main simulation script
Dependencies
Install required packages:

BASH
pip install numpy scipy matplotlib
File Descriptions
1. components.py
Defines physical components in the beamline:

PYTHON
class Target:       # Particle target with energy loss calculations
class Aperture:     # Beam-limiting aperture
class Magnets:      # Magnetic transport system
class Focalplane:   # Detection plane geometry
Key Features:

Energy loss calculations using interpolation from ESP.dat
Support for arbitrary focal plane geometries
Magnetic field coefficients loaded from external files
2. processes.py
Core simulation processes:

PYTHON
Beam_init()    # Generates particle beam with energy/angle distributions
Beam_trans()   # Transports particles through magnetic elements
Beam_hit()     # Tracks particle impacts on focal plane
Physics Models:

Random particle generation in target material
Energy loss via stopping power (ESP)
Magnetic transport using transfer matrices
Geometric intersection calculations
3. syn_MPR.py
Main simulation script:

PYTHON
# Initialize components
target = Target('CH2', r=0.02, thickness=0.00016)
aperture = Aperture('circular', r=0.022, distance=0.2)
Magnet = Magnets('TM.txt', reference_energy=14, length=0.8)

# Generate and transport beam
Beam = Beam_init('generate', target, aperture)
Beam = Beam_trans(Magnet, Beam)

# Detect impacts and visualize
focalplane = Focalplane('normal', position=0.86)
record = Beam_hit(Magnet, focalplane, Beam)
Execution
Run the simulation:

BASH
python syn_MPR.py
Expected Output:

Terminal logs of particle generation/transport
Interactive matplotlib window showing:
Particle impact positions on focal plane
Scatter plot of (l, y) coordinates
Example Output

Configuration
Modify these parameters in syn_MPR.py:

PYTHON
# Target properties
target = Target('CH2', r=0.02, thickness=0.00016)

# Beam energy range (prompted during execution)
Emin = ...  # Minimum energy [MeV]
Emax = ...  # Maximum energy [MeV]

# Magnet settings
Magnet = Magnets('TM.txt', reference_energy=14, length=0.8)

# Focal plane type
focalplane = Focalplane('normal', position=0.86)
# focalplane = Focalplane('arbitrary', ...)
Data Files
ESP.dat - Stopping power data (MeV/(g/cm²))
Format: Energy[MeV] Stopping_Power
TM.txt - Magnet transfer matrix coefficients
Format per line: C1 C2 C3 C4 C5 p1 p2 p3 p4 p5 p6
Where p1-p6 are exponent flags (0/1)
Key Features
🎯 Realistic beam generation with energy loss
🧲 Magnetic transport using transfer matrices
📐 Support for arbitrary focal plane geometries
📊 Interactive visualization of results
⚡ Optimized using NumPy and SciPy
For questions or contributions, contact [your email/contact info].
