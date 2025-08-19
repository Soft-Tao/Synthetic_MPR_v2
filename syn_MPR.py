import numpy as np
import matplotlib.pyplot as plt
from components import Target, Aperture, Magnets, Focalplane
from processes import Beam_init, Beam_trans, Beam_hit

target = Target(type='CH2', thickness=0.00016, shape="circle", geometry=[0.02])
#target = Target(type='CH2', thickness=0.00016, shape="rectangle", geometry=[0.03, 0.04189])

#aperture = Aperture(distance=0.2, shape="circle", geometry=[0.022])
aperture = Aperture(distance=0.2, shape="rectangle", geometry=[0.03, 0.04])

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

# histogram of l
plt.hist(x, bins=100)

# calculate FWHM （full width at half maximum） from histogram
hist, bin_edges = np.histogram(x, bins=100)
max_value = max(hist)
for i in range(len(hist)):
    if hist[i] >= max_value/2:
        break
for j in range(len(hist)):
    if hist[-j] >= max_value/2:
        break
FWHM = - bin_edges[i] + bin_edges[-j]

# calculate average from hist
sum = 0
for i in range(len(hist)):
    sum += hist[i] * bin_edges[i]
average = sum / len(x)
# calculate standard deviation from hist
sum = 0
for i in range(len(hist)):
    sum += hist[i] * (bin_edges[i] - average)**2
standard_deviation = np.sqrt(sum / len(x))

plt.xlabel('l [m]')
plt.ylabel('Counts')
plt.title(f'average = {average:.4f} m, FWHM = {FWHM:.4f} m, standard deviation = {standard_deviation:.4f} m')
plt.show()