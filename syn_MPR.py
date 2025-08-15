import numpy as np
import matplotlib.pyplot as plt
from components import Target, Aperture, Magnets, Focalplane
from processes import Beam_init, Beam_trans, Beam_hit

target = Target('CH2', r=0.02, thickness=0.00016)
aperture = Aperture('circular', r=0.022, distance=0.2)

Beam = Beam_init('generate', target, aperture)

'''# plot histgram of energy
import matplotlib.pyplot as plt
plt.hist([i[0] for i in Beam], bins=100)
plt.xlabel('Energy [MeV]')
plt.ylabel('Counts')
plt.show()'''

Magnet = Magnets('TM.txt', reference_energy=14, length=0.8)

Beam = Beam_trans(Magnet, Beam)

focalplane = Focalplane('normal', position=0.86)
#focalplane = Focalplane('arbitrary', position=0.86, geometry=[(-0.5, 0.81), (-0.324, 0.81), (-0.216, 0.84), (-0.072, 0.86), (0.042, 0.86), (0.126, 0.84), (0.5, 0.82)])

record = Beam_hit(Magnet, focalplane, Beam)

# scatter record
x = [i[0] for i in record]
y = [i[1] for i in record]
plt.scatter(x, y, s=2)
plt.xlabel('l [m]')
plt.ylabel('y [m]')
plt.axis('equal')
plt.show()