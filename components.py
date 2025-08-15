import numpy as np
from scipy.interpolate import interp1d

class Target:
    def __init__(self, type: str, thickness: float, r: float):
        self.type = type
        self.thickness = thickness
        self.r = r

        self.ESP = None
        self.density = None

        if self.type == 'CH2':
            self.density = 0.5 # [g/cm^3] range: 0.1-0.96

            with open("ESP.dat", 'r') as f:
                lines = f.readlines()
                ESP = np.array([line.strip().split() for line in lines], dtype=float)
            self.ESP = ESP # [MeV/(g/cm^2)]
            f.close()
        
        self.ESP_interp = interp1d(self.ESP[:,0], self.ESP[:,1])
    
    def get_ESP(self, E):
        return self.ESP_interp(E)
class Aperture:
    def __init__(self, type: str, r: float, distance: float):
        self.type = type
        self.distance = distance
        self.r = r


class Magnets:
    def __init__(self, file_path: str, reference_energy: float, length: float): # [MeV]
        self.reference_energy = reference_energy
        self.TM = []
        self.length = length

        with open(file_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                TMi = []
                for i in line.strip().split()[:5]:
                    TMi.append(float(i))
                for l in line.strip().split()[5]:
                    TMi.append(int(l))
                self.TM.append(TMi)

class Focalplane:
    def __init__(self, type: str, position: float, geometry = None): # geometry: (x, z) list of tuples, z = 0 @ position
        self.type = type
        self.position = position

        if self.type == 'normal':
            self.geometry = [(x, self.position, x+0.5) for x in np.linspace(-0.5, 0.5, 101)]
        elif self.type == 'arbitrary':
            if geometry == None:
                raise Exception("you CANNOT generate an arbitrary focal plane with no geometry specified!")
            else:
                l = 0
                self.geometry = [(geometry[0][0], geometry[0][1], 0)]
                for i in range(1, len(geometry)):
                    l += np.sqrt((geometry[i][0] - geometry[i-1][0])**2 + (geometry[i][1] - geometry[i-1][1])**2)
                    self.geometry.append((geometry[i][0], geometry[i][1], l))