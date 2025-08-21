import numpy as np
import matplotlib.pyplot as plt
from components import Target, Aperture, Magnets, Focalplane
from processes import Beam_init, Beam_trans, Beam_hit

#target = Target(type='CH2', thickness=0.00016, shape="circle", geometry=[0.02])
#target = Target(type='CH2', thickness=0.00016, shape="rectangle", geometry=[0.0271, 0.0579])
target = Target(type='CH2', thickness=0.00016, shape="rectangle", geometry=[0.0194, 0.0647])

aperture = Aperture(distance=0.2, shape="circle", geometry=[0.022])
#aperture = Aperture(distance=0.2, shape="rectangle", geometry=[0.0214, 0.0712])

Magnet = Magnets('TM.txt', reference_energy=14, length=1.75)

focalplane = Focalplane('normal', position=1.81)
#geo_lst = [(x, 1.67 + 1.8 * x) for x in np.linspace(-0.2, 0.3, 21)]

#focalplane = Focalplane('arbitrary', position=0.86, geometry=geo_lst)

Ave_lst = []
Fwhm_lst = []
Energy_lst = np.linspace(9.5, 19.5, 21)

for ene in Energy_lst:

    Beam = Beam_init('generate', target, aperture, Emin = ene, Emax = ene, Nbins = 1, Npart = 7000)
    Beam = Beam_trans(Magnet, Beam)
    record = Beam_hit(Magnet, focalplane, Beam)

    # scatter record
    x = [i[0] for i in record]
    y = [i[1] for i in record]
    '''plt.scatter(x, y, s=2)
    plt.xlabel('l [m]')
    plt.ylabel('y [m]')
    plt.axis('equal')
    plt.show()'''

    # histogram of l
    #plt.hist(x, bins=np.linspace(0, 1, 21))

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
    
    Ave_lst.append(average)
    Fwhm_lst.append(FWHM)

'''plt.xlabel('l [m]')
plt.ylabel('Counts')
plt.title(f'average = {average:.4f} m, FWHM = {FWHM:.4f} m, standard deviation = {standard_deviation:.4f} m')
plt.show()'''

plt.plot(Energy_lst, Ave_lst)
plt.xlabel('Energy [MeV]')
plt.ylabel('average [m]')
plt.show()

plt.scatter(Energy_lst, Fwhm_lst, c = 'black')
plt.plot([min(Energy_lst), max(Energy_lst)], [0.05, 0.05], 'k--')
plt.xlabel('Energy [MeV]')
plt.ylabel('FWHM [m]')
plt.show()

ene_resolution = []
for i in range(len(Energy_lst)):
    if (i-1) % 2 == 0:
        ene_resolution.append([Energy_lst[i], (Fwhm_lst[i]*(Energy_lst[i+1] - Energy_lst[i-1])/(Ave_lst[i+1] - Ave_lst[i-1]))/Energy_lst[i]*100])
# plot ene_resolution
plt.scatter([i[0] for i in ene_resolution], [i[1] for i in ene_resolution], c = 'black')
plt.xlabel('Energy [MeV]')
# plot y = 4
plt.plot([min(Energy_lst), max(Energy_lst)], [4, 4], 'k--')
plt.ylabel(r'Energy resolution [%]')
plt.show()
    