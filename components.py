import numpy as np
from scipy.interpolate import interp1d, LinearNDInterpolator

def load_numeric_table(fname):
    start = None
    with open(fname, "r") as f:
        for i, line in enumerate(f):
            parts = line.strip().split()
            if not parts:
                continue
            try:
                _ = [float(x) for x in parts]
                start = i
                break
            except ValueError:
                continue
    return np.loadtxt(fname, skiprows=start)

def sample_u_vectorized(f, N):
    """
    Generate N*u in [0, 1] with probability density ∝ u * f(2u^2 - 1)
    """
    M = f(1.0)
    u_samples = []

    while len(u_samples) < N:
        u_candidate = np.random.uniform(0, 1, size=5 * N)
        w = u_candidate * f(2 * u_candidate ** 2 - 1)
        accept = u_candidate[np.random.uniform(0, M, size=u_candidate.size) < w]
        u_samples.extend(accept.tolist())
    return np.array(u_samples[:N])

def sample_u_vectorized_with_umin(f, g, N, u_min):
    """
    Generate N*u in [u_min, 1] with probability density ∝ u * f(2u^2 - 1)

    Parameters:
    f: probability density function
    g: integral of f from 0 to x
    N: number of samples needed
    u_min: minimum value of u, u ∈ [u_min, 1]

    return:
    u_samples: sampled u array
    prob_ratio: probability ratio of [u_min,1] to [0,1]
    """
    M = f(1.0)
    u_samples = []

    # calculate prob_ratio
    gp = lambda uu: 1/4.0 * (g(2*uu**2-1) - g(-1))
    prob_ratio = (gp(1.0) - gp(u_min)) / gp(1.0)

    while len(u_samples) < N:
        u_candidate = np.random.uniform(u_min, 1.0, size=5 * N)
        w = u_candidate * f(2 * u_candidate ** 2 - 1)
        accept = u_candidate[np.random.uniform(0, M, size=u_candidate.size) < w]
        u_samples.extend(accept.tolist())

    return np.array(u_samples[:N]), prob_ratio


class Target:
    def __init__(self, type: str, mass_thickness: float, shape: str, geometry: list=None, area=None, H_W_ratio=None,
                 density=0.94, ESP_file='ESP(PSTAR).dat'):
        # geometry: [r] if shape == "circle", [w, h] if shape == "rectangle"
        """
        Initialize a Target object.

        Parameters:
        type (str): Material type of the target, e.g., "CH2".
        mass_thickness (float): mass density of the target, unit: [mg/cm^2]
        thickness (float): Thickness of the target. [m]
        shape (str): Shape of the target, e.g., "circle", "rectangle".
        geometry (list): List describing geometric dimensions of the target.
                         The meaning depends on the shape:
                           - If shape = "disk", geometry = [radius]
                           - If shape = "square", geometry = [side length]
                           - If shape = "rectangle", geometry = [width, height]
        """
        self.type = type
        self.mass_thickness = mass_thickness
        self.shape = shape
        if geometry is not None:
            self.geometry = geometry
            self._get_area()
        elif area is not None and H_W_ratio is not None:
            self.area = area
            self.H_W_ratio = H_W_ratio
            self._get_geometry()

        self.ESP = None
        self.density = None

        if self.type == 'CH2':
            self.density = density  # [g/cm^3] range: 0.1-0.96
            ESP = load_numeric_table(ESP_file)
            self.ESP = ESP
        else:
            raise Exception("Target: Invalid type!")

        self.thickness = self.mass_thickness / self.density * 1E-5  # [m]
        self.ESP_interp = interp1d(self.ESP[:, 0], self.ESP[:, 1])
        self._get_max_r()

    def _get_area(self):
        if self.shape == 'circle':
            self.area = np.pi * self.geometry[0]**2
            self.H_W_ratio = 1.0
        elif self.shape == 'rectangle':
            self.area = self.geometry[0] * self.geometry[1]
            self.H_W_ratio = self.geometry[1] / self.geometry[0]

    def _get_geometry(self):
        if self.shape == 'circle':
            r = np.sqrt(self.area / np.pi)
            self.geometry = [r]
        elif self.shape == 'rectangle':
            width = np.sqrt(self.area / self.H_W_ratio)
            height = width * self.H_W_ratio
            self.geometry = [width, height]

    def _get_max_r(self):
        if self.shape == 'circle':
            self.max_r = self.geometry[0]
        elif self.shape == 'rectangle':
            self.max_r = 0.5 * np.sqrt(self.geometry[0]**2 + self.geometry[1]**2)

    def get_ESP(self, E):
        return self.ESP_interp(E)

    def get_initialPosition(self):
        if self.shape == "circle":
            while (True):
                x = np.random.uniform(-self.geometry[0], self.geometry[0])
                y = np.random.uniform(-self.geometry[0], self.geometry[0])
                if x ** 2 + y ** 2 <= self.geometry[0] ** 2:
                    break
            z = np.random.uniform(0, self.thickness)
            return x, y, z
        elif self.shape == "rectangle":
            return np.random.uniform(-self.geometry[0] / 2, self.geometry[0] / 2), np.random.uniform(
                -self.geometry[1] / 2, self.geometry[1] / 2), np.random.uniform(0, self.thickness)
        else:
            raise Exception("Target: Invalid shape!")

    def get_initialPosition_N(self, N):
        if self.shape == "circle":
            x_list, y_list = [], []
            while len(x_list) < N:
                x = np.random.uniform(-self.geometry[0], self.geometry[0], N)
                y = np.random.uniform(-self.geometry[0], self.geometry[0], N)
                idx = x ** 2 + y ** 2 <= self.geometry[0] ** 2
                x_list.extend(x[idx])
                y_list.extend(y[idx])
            return np.array(x_list[0:N]), np.array(y_list[0:N])
        elif self.shape == "rectangle":
            return np.random.uniform(-self.geometry[0] / 2, self.geometry[0] / 2, N), np.random.uniform(
                -self.geometry[1] / 2, self.geometry[1] / 2, N)
        else:
            raise Exception("Target: Invalid shape!")

    def generate_particles(self, Np, f, g, u_min):
        if self.shape == "circle":
            x_list, y_list = [], []
            while len(x_list) < Np:
                batch_size = 2 * (Np - len(x_list))
                x_c = np.random.uniform(-self.geometry[0], self.geometry[0], batch_size)
                y_c = np.random.uniform(-self.geometry[0], self.geometry[0], batch_size)
                mask = x_c ** 2 + y_c ** 2 <= self.geometry[0] ** 2
                x_list.extend(x_c[mask].tolist())
                y_list.extend(y_c[mask].tolist())
            x_arr = np.array(x_list[:Np])
            y_arr = np.array(y_list[:Np])
        elif self.shape == "rectangle":
            x_arr = np.random.uniform(-self.geometry[0] / 2, self.geometry[0] / 2, Np)
            y_arr = np.random.uniform(-self.geometry[1] / 2, self.geometry[1] / 2, Np)
        else:
            raise NotImplementedError("Target: Invalid shape!")

        z_arr = np.random.uniform(0, self.thickness, Np)
        phi_arr = np.random.uniform(0, 2 * np.pi, Np)
        # from processes import sample_u_vectorized

        # u_arr = sample_u_vectorized(f, Np)
        u_arr, frac = sample_u_vectorized_with_umin(f, g, Np, u_min)

        return np.stack([x_arr, y_arr, z_arr, u_arr, phi_arr], axis=1), frac


class Aperture:
    def __init__(self, distance: float, shape: str, solid_angle=None, H_W_ratio=None, geometry: list=None):
        """
        Initialize an Aperture object.

        Parameters:
        distance: Distance between the Target and the Aperture [m]
        shape: 'circle', 'rectangle'
        solid_angle: solid angle [msr]
        H_W_ratio: H/W of the aperture
        geometry: [r] for 'circle', [w, h] for 'rectangle', [m]
        """
        self.distance = distance
        self.shape = shape
        if geometry is not None:
            self.geometry = geometry
            self._get_solid_angle()
        elif solid_angle is not None and H_W_ratio is not None:
            self.solid_angle = solid_angle
            self.H_W_ratio = H_W_ratio
            self._get_geometry()
        else:
            raise ValueError("ERROR occurred when defining aperture!")
        self._get_max_r()

    def _get_solid_angle(self):
        if self.shape == "circle":
            self.area = np.pi * self.geometry[0] ** 2
            self.solid_angle = self.area / self.distance**2
            self.H_W_ratio = 1.0
        elif self.shape == "rectangle":
            self.area = self.geometry[0] * self.geometry[1]
            self.solid_angle = self.area / self.distance**2
            self.H_W_ratio = self.geometry[1] / self.geometry[0]

    def _get_geometry(self):
        self.area = self.solid_angle * 1E-3 * self.distance**2
        if self.shape == 'circle':
            r = np.sqrt(self.area / np.pi)
            self.geometry = [r]
        elif self.shape == 'rectangle':
            width = np.sqrt(self.area / self.H_W_ratio)
            height = self.H_W_ratio * width
            self.geometry = [width, height]

    def _get_max_r(self):
        if self.shape == 'circle':
            self.max_r = self.geometry[0]
        elif self.shape == 'rectangle':
            self.max_r = 0.5 * np.sqrt(self.geometry[0]**2 + self.geometry[1]**2)

    def isPassed(self, x, y):
        if self.shape == "circle":
            return x ** 2 + y ** 2 <= self.geometry[0] ** 2
        elif self.shape == "rectangle":
            return np.logical_and(np.abs(x) <= self.geometry[0] / 2, np.abs(y) <= self.geometry[1] / 2)

class Magnets:
    def __init__(self, type: str, file_path: str, reference_energy: float = None):
        '''
        ===If type is 2d===
        Initialize a Magnets object using transfer map from COSY INFINITY.

        Parameters:
        file_path: path to the TM.txt file
        reference_energy: reference energy [MeV]

        ===If type is 3d===
        Initialize a Magnets object using three dimensional B-field data (e.g. from COMSOL).

        Parameters:
        file_path: path to the B-field data
        '''
        self.type = type
        if self.type == '2d':
            if reference_energy is None:
                raise Exception("you must provide [reference_energy] if you're trying to initiate a 2d magnets!")
            self.reference_energy = reference_energy
            self.TM = []

            with open(file_path, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    TMi = []
                    for i in line.strip().split()[:5]:
                        TMi.append(float(i))
                    for l in line.strip().split()[5]:
                        TMi.append(int(l))
                    self.TM.append(TMi)
            self.TM = np.array(self.TM)
        if self.type == '3d':
            with open(file_path, 'r') as f:
                lines = f.readlines()
                data = []
                for line in lines:
                    if line[0] != "%":
                        data.append([float(i) for i in line.strip().split(',')])
            f.close()
            data = np.array(data)
            self.xgrid = data[:, 0]
            self.ygrid = data[:, 1]
            self.zgrid = data[:, 2]
            self.Bx = data[:, 3]
            self.By = data[:, 4]
            self.Bz = data[:, 5]
            self._Bx_interp = None
            self._By_interp = None
            self._Bz_interp = None
        
    def _initialize_interpolators(self):
        self._Bx_interp = LinearNDInterpolator((self.xgrid, self.ygrid, self.zgrid), self.Bx)
        self._By_interp = LinearNDInterpolator((self.xgrid, self.ygrid, self.zgrid), self.By)
        self._Bz_interp = LinearNDInterpolator((self.xgrid, self.ygrid, self.zgrid), self.Bz)
    
    def get_B_at(self, X):
        if self._Bx_interp is None:
            self._initialize_interpolators()
        Bx = np.array(self._Bx_interp(X))
        By = np.array(self._By_interp(X))
        Bz = np.array(self._Bz_interp(X))
        return np.stack((Bx, By, Bz), axis=1)

class Focalplane:
    def __init__(self, type: str, position: float=None, geometry: list=None, length=1):
        '''
        Initialize a Focalplane object.

        Parameters:
        **In 2d case (COSY TM based)**
        type: 'normal' or 'arbitrary'
        position: focal plane position [m] behind magnets' exit
        geometry: list (x, z) of tuples, z = 0 @ magnets' exit

        **In 3d case (e.g. COMSOL based)**
        type: 'arbitrary' only
        geometry: list (x, y) of tuples, origin @ target's geometric center
        (position parameter is ignored in 3d case)
        '''
        self.type = type
        self.position = position

        if self.type == 'normal':
            if self.position is None:
                raise Exception("you must provide [position] if you're trying to initiate a normal fp!")
            self.geometry = [(x, self.position, x + length/2.) for x in np.linspace(-length/2., length/2., 11)]
            self.geometry = np.array(self.geometry)
        elif self.type == 'arbitrary':
            if geometry is None:
                raise Exception("you CANNOT generate an arbitrary focal plane with no geometry specified!")
            else:
                l = 0
                self.geometry = [(geometry[0][0], geometry[0][1], 0)]
                for i in range(1, len(geometry)):
                    l += np.sqrt(
                        (geometry[i][0] - geometry[i - 1][0]) ** 2 + (geometry[i][1] - geometry[i - 1][1]) ** 2)
                    self.geometry.append((geometry[i][0], geometry[i][1], l))
                self.geometry = np.array(self.geometry)
