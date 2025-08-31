import numpy as np
from scipy.interpolate import interp1d



def load_numeric_table(fname):
    start = None
    with open(fname, "r") as f:
        for i, line in enumerate(f):
            parts = line.strip().split()
            if not parts:
                continue
            try:
                # 全部都能转为 float，才算有效数据行
                _ = [float(x) for x in parts]
                start = i
                break
            except ValueError:
                continue
    return np.loadtxt(fname, skiprows=start)


def sample_u_vectorized(f, N):
    """
    生成 N 个 u ∈ [0,1]，概率密度 ∝ u * f(2u^2 - 1)
    向量化批量生成
    """
    M = f(1.0)  # 最大值
    u_samples = []

    while len(u_samples) < N:
        # 生成候选 u
        u_candidate = np.random.uniform(0, 1, size=5 * N)
        w = u_candidate * f(2 * u_candidate ** 2 - 1)
        accept = u_candidate[np.random.uniform(0, M, size=u_candidate.size) < w]
        u_samples.extend(accept.tolist())
    return np.array(u_samples[:N])


def sample_u_vectorized_with_umin(f, g, N, u_min):
    """
    生成 N 个 u ∈ [u_min, 1]，概率密度 ∝ u * f(2u^2 - 1)

    :param f: 概率密度函数
    :param g: 代表 f 从 0 到 x 的积分函数
    :param N: 需要采样的数量
    :param u_min: 最小取值，u ∈ [u_min, 1]

    :return: u_samples: 采样得到的 u 数组
            prob_ratio: [u_min,1]区间概率与[0,1]区间概率的比值
    """
    M = f(1.0)  # 最大值
    u_samples = []

    # 计算 [u_min, 1] 区间概率与 [0, 1] 区间概率的比值
    gp = lambda uu: 1/4.0 * (g(2*uu**2-1) - g(-1))
    prob_ratio = (gp(1.0) - gp(u_min)) / gp(1.0)

    while len(u_samples) < N:
        # 生成候选 u，范围在 [u_min, 1]，通过线性变换映射
        u_candidate = np.random.uniform(u_min, 1.0, size=5 * N)
        w = u_candidate * f(2 * u_candidate ** 2 - 1)

        # 接受的标准
        accept = u_candidate[np.random.uniform(0, M, size=u_candidate.size) < w]

        # 扩展样本列表
        u_samples.extend(accept.tolist())

    # 取前 N 个采样值
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
                           - If shape = "rectangle", geometry = [radius, height]
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
            # with open("ESP.dat", 'r') as f:
            #     lines = f.readlines()
            #     ESP = np.array([line.strip().split() for line in lines], dtype=float)
            ESP = load_numeric_table(ESP_file)
            self.ESP = ESP  # [MeV/(g/cm^2)]
            # f.close()
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
            x_arr = np.random.uniform(-self.geometry[0] / 2, self.geometry[0] / 2)
            y_arr = np.random.uniform(-self.geometry[1] / 2, self.geometry[1] / 2)
        else:
            raise NotImplementedError("Target: Invalid shape!")

        z_arr = np.random.uniform(0, self.thickness, Np)
        phi_arr = np.random.uniform(0, 2 * np.pi, Np)
        # from processes import sample_u_vectorized

        # u_arr = sample_u_vectorized(f, Np)
        u_arr, frac = sample_u_vectorized_with_umin(f, g, Np, u_min)

        return np.stack([x_arr, y_arr, z_arr, u_arr, phi_arr], axis=1), frac


class Aperture:
    def __init__(self, distance: float, shape: str, solid_angle=None, H_W_ratio=None, geometry=None):
        """
        :param distance: Distance between the Target and the Aperture [m]
        :param shape: 'circle', 'rectangle'
        :param solid_angle: solid angle [msr]
        :param H_W_ratio: H/W of the aperture
        :param geometry: [r] for 'circle', [w, h] for 'rectangle', [m]
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
            return (np.abs(x) <= self.geommetry[0] / 2 and np.abs(y) <= self.geometry[1] / 2)


class Magnets:
    def __init__(self, file_path: str, reference_energy: float):    # , length: float):  # [MeV]
        self.reference_energy = reference_energy
        self.TM = []
        # self.length = length

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


class Focalplane:
    def __init__(self, type: str, position: float, geometry=None, length=1):  # geometry: (x, z) list of tuples, z = 0 @ position
        self.type = type
        self.position = position

        if self.type == 'normal':
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
