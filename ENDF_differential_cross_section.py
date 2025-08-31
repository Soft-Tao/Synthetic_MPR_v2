import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre
import endf
from scipy.interpolate import interp1d

def read_ang_dist_data(filename):
    mat = endf.Material(filename)
    E = np.array(mat.section_data[4, 2]['legendre']['E'])
    a_l = np.array(mat.section_data[4, 2]['legendre']['a_l'])
    ones_col = np.ones((a_l.shape[0], 1))
    # 水平拼接：在每一行前加一个 1
    a_l_new = np.hstack((ones_col, a_l))
    return E, a_l_new
    # return E, a_l

def read_sigma_data(filename):
    mat = endf.Material(filename)
    sigma = mat.section_data[3, 2]['sigma']
    return sigma.x, sigma.y

def find_closest_index(array, value):
    idx = np.abs(array - value).argmin()
    return idx

def compute_f(theta, coeffs):
    # theta是弧度数组，coeffs是a_l列表
    cos_theta = np.cos(theta)
    f = np.zeros_like(theta)
    for l, a_l in enumerate(coeffs):
        P_l = legendre(l)  # 生成l阶勒让德多项式
        print(f"l = {l}, a_l = {a_l:.2E}, coeff = {(2*l + 1)/2 * a_l:.2e}, P_l(0) = {P_l(1)}")
        f += (2*l + 1)/2 * a_l * P_l(cos_theta)
    return f
def make_f_v(coeffs):
    """
    输入: coeffs 系数列表 [a_0, a_1, ..., a_l]
    输出: f(v) 函数, v=cos(theta)
    """
    coeffs = np.asarray(coeffs)
    def f(v):
        v_arr = np.atleast_1d(v)
        result = np.zeros_like(v_arr)
        for l, a_l in enumerate(coeffs):
            P_l = np.polynomial.legendre.Legendre.basis(l)
            result += (2*l + 1)/2 * a_l * P_l(v_arr)
        return float(result[0]) if np.isscalar(v) else result
    return f

def int_P(l):
    if l > 0:
        return (-np.polynomial.legendre.Legendre.basis(l-1) + np.polynomial.legendre.Legendre.basis(l+1) +
            np.polynomial.legendre.Legendre.basis(l-1)(0) - np.polynomial.legendre.Legendre.basis(l+1)(0)) / (2*l + 1.0)
    else:
        return lambda x:x

def make_g_v(coeffs):
    """
    输入: coeffs 系数列表 [a_0, a_1, ..., a_l]
    输出: f(v) 函数, v=cos(theta)
    """
    coeffs = np.asarray(coeffs)
    def g(v):
        v_arr = np.atleast_1d(v)
        result = np.zeros_like(v_arr, dtype=float)
        for l, a_l in enumerate(coeffs):
            int_P_l = int_P(l)
            result += (2*l + 1)/2 * a_l * int_P_l(v_arr)
        return float(result[0]) if np.isscalar(v) else result
    return g


class differential_cross_section:
    def __init__(self, filename):
        energies, coeffs_list = read_ang_dist_data(filename)
        E, sigma_s = read_sigma_data(filename)
        self.E_list = energies
        self.coeffs_list = coeffs_list
        self.sigma_s_list = sigma_s
        self._interp_sigma_s = None
        self._interps_coeffs = None
    def compute_sigma_s(self, E):
        # E[eV] -> sigma_s[barn] or [1E-28 m^2]
        if self._interp_sigma_s is None:
            self._interp_sigma_s = interp1d(
                self.E_list, self.sigma_s_list,
                kind='linear', fill_value="extrapolate"
            )
        return float(self._interp_sigma_s(E))

    def compute_coeffs(self, E):
        if self._interps_coeffs is None:
            self._interps_coeffs = [
                interp1d(self.E_list, self.coeffs_list[:, i],
                         kind='linear', fill_value="extrapolate")
                for i in range(self.coeffs_list.shape[1])
            ]
        return np.array([interp(E) for interp in self._interps_coeffs])

    def compute_f(self, E):
        coeffs = self.compute_coeffs(E)
        return make_f_v(coeffs)

    def compute_g(self, E):
        coeffs = self.compute_coeffs(E)
        return make_g_v(coeffs)

if __name__ == "__main__":
    filename = 'E4R84432_e4.endf.endf2gnd.endf'
    energy_input = 14E6     # eV
    cross_section = differential_cross_section(filename)

    plt.figure(figsize=(5,4))
    plt.plot(cross_section.E_list, cross_section.sigma_s_list)  # x轴用角度（度）表示
    plt.ylabel('$\\sigma_{\\text{s}}\\ [\\text{barns}]$')
    plt.xlabel('$E_{\\text{n}}\\ [\\text{eV}]$')
    plt.title('n-p Scattering $\\sigma_{\\text{s}}$')
    # plt.xlim((0,90))
    # plt.ylim(bottom=0)
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    # plt.savefig('all.png',dpi=300)
    # plt.savefig('all.pdf',dpi=300)
    plt.show()
    #
    # idx_14 = np.argmin(np.abs(cross_section.E_list - 14.0E6))
    # print(f"{cross_section.sigma_s_list[idx_14]} barns")



    f_v = cross_section.compute_f(energy_input)
    sigma_s = cross_section.compute_sigma_s(energy_input)
    v_list = np.linspace(-1, 1, 100)
    beta_list = np.arccos(v_list)
    theta_p_lab_list = np.pi / 2 - beta_list / 2
    theta_n_lab_list = beta_list / 2
    f_list = f_v(v_list)
    f_p_lab = sigma_s / (2.0*np.pi) * f_list * 4 * np.cos(theta_p_lab_list)
    f_n_lab = sigma_s / (2.0*np.pi) * f_list * 4 * np.cos(theta_n_lab_list)

    plt.figure(figsize=(5,4))
    plt.plot(np.rad2deg(theta_p_lab_list), f_p_lab, label ='recoil proton')  # x轴用角度（度）表示
    plt.plot(np.rad2deg(theta_n_lab_list), f_n_lab, label ='scattered neutron')  # x轴用角度（度）表示
    plt.ylabel('$\\text{d}\\sigma/\\text{d}\\Omega_{\\text{LAB}}\\ [\\text{barns}/\\text{sr}]$')
    plt.xlabel('$\\theta_{\\text{p}}\\ [\\text{degrees}]$')
    plt.title('n-p Scattering $\\text{d}\\sigma/\\text{d}\\Omega_{\\text{LAB}},\\ \\ E_{\\text{n}}\\,=\\,$'+f'{energy_input/1E6:.1f} MeV')
    plt.xlim((0,90))
    plt.ylim(bottom=0)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    # plt.savefig('all.png',dpi=300)
    # plt.savefig('all.pdf',dpi=300)
    plt.show()

    u_n_lab_list = np.cos(theta_n_lab_list)
    u_p_lab_list = np.cos(theta_p_lab_list)
    plt.figure(figsize=(5,4))
    plt.plot(u_p_lab_list, f_p_lab, label ='recoil proton')
    plt.plot(u_n_lab_list, f_n_lab, label ='scattered neutron')
    plt.ylabel('$\\text{d}\\sigma/\\text{d}\\Omega_{\\text{LAB}}\\ [\\text{barns}/\\text{sr}]$')
    plt.xlabel('$\\cos\\theta_{\\text{p}}$')
    plt.title('n-p Scattering $\\text{d}\\sigma/\\text{d}\\Omega_{\\text{LAB}},\\ \\ E_{\\text{n}}\\,=\\,$'+f'{energy_input/1E6:.1f} MeV')
    plt.xlim((0,1))
    plt.ylim(bottom=0)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    # plt.savefig('all.png',dpi=300)
    # plt.savefig('all.pdf',dpi=300)
    plt.show()
