import numpy as np
import qutip
from .constants import *
from scipy.linalg import eig

def adiabatic_passage_f_sweep(atom, F, f_res, f_range, B_0, B_rf):
    B_homogenous = B_0*1e-4
    B_int = B_rf*1e-4

    g_F = atom.g_J*(F*(F+1) - atom.I*(atom.I+1) + atom.J*(atom.J+1)) / (2*F*(F+1)) +\
        atom.g_I*(F*(F+1) + atom.I*(atom.I+1) - atom.J*(atom.J+1)) / (2*F*(F+1))

    omega_l_2 = g_F * mu_bohr * B_homogenous / hbar
    omega_q_2 = ((atom.g_J-atom.g_I)**2 * mu_bohr**2 * B_homogenous**2/(4*hbar* h*2*np.pi*atom.delta_E_hf)) 
    eigenvalues = np.zeros((f_res, int(2*F+1)))

    delta_range = np.linspace(-f_range/2, f_range/2, f_res)
    omega_rf_range = omega_l_2 + delta_range

    omega_rabi = mu_bohr * g_F * B_int / hbar
    H_int = hbar*omega_rabi/2 * np.mat(qutip.jmat(F, "x"))

    for i, delta in enumerate(delta_range):
        H_ze_l = hbar * delta * np.mat(qutip.jmat(F, "z"))
        H_ze_q = hbar * omega_q_2 * np.mat(1-(2*qutip.jmat(F, "z")/(2*atom.I+1))**2)
        H = H_ze_l + H_ze_q + H_int
        eigenvalues[i] = sorted(np.real(eig(H, right=False)))
    return eigenvalues, omega_rf_range/(2*np.pi)

def adiabatic_passage_B_sweep(atom, F, B_res, B_range, f, B_rf):
    B_int = B_rf * 1e-4

    g_F = atom.g_J*(F*(F+1) - atom.I*(atom.I+1) + atom.J*(atom.J+1)) / (2*F*(F+1)) +\
        atom.g_I*(F*(F+1) + atom.I*(atom.I+1) - atom.J*(atom.J+1)) / (2*F*(F+1))

    B_0_center = 2*np.pi*f*hbar/(g_F*mu_bohr) 

    omega_rabi = mu_bohr * g_F * B_int / hbar
    H_int = omega_rabi/2 * np.mat(qutip.jmat(F, "x"))

    eigenvalues = np.zeros((B_res, int(2*F+1)))

    B_array = B_0_center + np.linspace(-B_range/2, B_range/2, B_res)*1e-4

    for i, B_0 in enumerate(B_array):
        omega_l_2 = g_F * mu_bohr * B_0 / hbar
        omega_q_2 = ((atom.g_J-atom.g_I)**2 * mu_bohr**2 * B_0**2/(4*hbar* h*2*np.pi*atom.delta_E_hf)) 

        delta_f = omega_l_2 - f*2*np.pi

        H = delta_f * np.mat(qutip.jmat(F, "z"))
        H += omega_q_2 * np.mat(1-(2*qutip.jmat(F, "z")/(2*atom.I+1))**2)
        H += H_int

        eigenvalues[i] = sorted(np.real(eig(H)[0]))
    return eigenvalues, B_array
