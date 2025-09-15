import numpy as np
import matplotlib.pyplot as plt
from components import Target, Aperture, Magnets, Focalplane
from processes import Beam
import ENDF_differential_cross_section

class MPR:
    '''
    An integrated MPR system initialized by specifying the 4 components.
    It has following functions:
    1. calculate l-E relation, energy resolution in a given energy range.
    2. find optimal straight focal plane position (position and tilt angle).
    3. calculate response matrix of the system.
    '''
    def __init__(self, target: Target, aperture: Aperture, magnet: Magnets, focalplane: Focalplane):
        self.target = target
        self.aperture = aperture
        self.magnet = magnet
        self.focalplane = focalplane

    def performance(self, E_range = np.linspace(9, 21, 13), N_part = 50000, plot_save = False):
        output = []
        cross_section = ENDF_differential_cross_section.differential_cross_section(
            'E4R84432_e4.endf.endf2gnd.endf')

        for i, ene in enumerate(E_range):
            print(f"[{i}/{len(E_range)}] energy: {ene:.1f} MeV")
            beam = Beam.generate(self.target, self.aperture, ene, N_part, cross_section=cross_section)
            beam_transported = beam.trans(self.magnet)
            record = beam_transported.hit(self.focalplane)
            output.append([ene, record.l_mean, record.std_dev_l])
        
        output = np.array(output)
        plt.figure()
        plt.plot(output[:, 0], output[:, 1], '-s', color='k')
        plt.scatter(output[:, 0], output[:, 1], s=2, color='k')
        plt.title('L-E relation')
        plt.xlabel('energy [MeV]')
        plt.ylabel(r'$\bar{l}$ [m]')
        if plot_save:
            plt.savefig('L_E_relation.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        resolution = output[1:-1, 2]/(output[2:, 1] - output[:-2, 1])*(output[2:, 0] - output[:-2, 0])/E_range[1:-1]*100
        plt.figure()
        plt.plot(E_range[1:-1], resolution, '-s', color='k')
        plt.scatter(E_range[1:-1], resolution, s=2, color='k')
        plt.title('Energy resolution') 
        plt.xlabel('energy [MeV]')
        plt.ylabel(r'Energy resolution [%]')
        plt.ylim(ymin=0)
        if plot_save:
            plt.savefig('energy_resolution.png', dpi=300, bbox_inches='tight')
        plt.show()


target = Target(type='CH2', mass_thickness=15, shape="circle", area=0.001, H_W_ratio=1.0)

aperture = Aperture(distance=0.2, shape="circle", solid_angle=40, H_W_ratio=1.0)

Magnet = Magnets('beam_trans_benchmark/TM.txt', reference_energy=14)

focalplane = Focalplane('normal', position=0.256)

mpr = MPR(target, aperture, Magnet, focalplane)

mpr.performance()