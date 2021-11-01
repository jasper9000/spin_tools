import numpy as np
from collections.abc import Iterable
import qutip
from .constants import *


def E_hf_breit_rabi_legacy(F, mf, B, atom):
    x = (mu_bohr*(atom.g_J - atom.g_I)*B) / (h*atom.delta_E_hf)
    if F == 7/2:
        # F = I - S
        E_hf_tmp = -h*atom.delta_E_hf/(2*(2*atom.I+1))
        E_hf_tmp += atom.g_I*mu_bohr*mf*B
        E_hf_tmp -= h*atom.delta_E_hf/2 * \
            np.sqrt((1 + 4*mf*x/(2*atom.I+1) + x**2))
    else:
        # F = I + S
        E_hf_tmp = -h*atom.delta_E_hf/(2*(2*atom.I+1))
        E_hf_tmp += atom.g_I*mu_bohr*mf*B
        if np.abs(mf) == (atom.I+0.5):
            E_hf_tmp += h*atom.delta_E_hf/2 * (1 + np.sign(mf)*x)
        else:
            E_hf_tmp += h*atom.delta_E_hf/2 * \
                np.sqrt((1 + 4*mf*x/(2*atom.I+1) + x**2))
    return E_hf_tmp
    
def E_hf_breit_rabi(B, atom, F):
    # convert to Tesla
    B_tmp = B * 1e-4

    if (atom.J != 1/2):
        raise(ValueError("The Breit-Rabi formula can only be applied to a manifold with J=1/2."))
    
    if (F != atom.I+atom.J) and (F != atom.I-atom.J):
        raise(ValueError("F != I +- J."))

    # make B iterable if it isnt already
    if not isinstance(B_tmp, Iterable):
        B_tmp = [B_tmp]

    mFs = np.arange(-F, F+1, 1, dtype=np.float)
    E_hf_tmp = np.zeros((len(mFs), len(B_tmp)))
    # E_hf_tmp = zeros_like(mFs, len(B))

    for i, b in enumerate(B_tmp):
        x = (mu_bohr*(atom.g_J - atom.g_I)*b / (h*atom.delta_E_hf))
        if (F == atom.I - 0.5):
            # F = I - J
            E_hf_tmp[:,i] += -h*atom.delta_E_hf/(2*(2*atom.I+1))
            E_hf_tmp[:,i] += atom.g_I*mu_bohr*mFs*b
            E_hf_tmp[:,i] -= h*atom.delta_E_hf/2 * \
                np.sqrt((1 + 4*mFs*x/(2*atom.I+1) + x**2))
        else:
            # F = I + J
            E_hf_tmp[:,i] += -h*atom.delta_E_hf/(2*(2*atom.I+1))
            E_hf_tmp[:,i] += atom.g_I*mu_bohr*mFs*b

            # cases
            E_hf_tmp[1:-1, i] += h*atom.delta_E_hf/2 * \
                    np.sqrt((1 + 4*mFs[1:-1]*x/(2*atom.I+1) + x**2))
            
            E_hf_tmp[0, i] += h*atom.delta_E_hf/2 * (1 + np.sign(mFs[0])*x)
            E_hf_tmp[-1, i] += h*atom.delta_E_hf/2 * (1 + np.sign(mFs[-1])*x)
    return E_hf_tmp

def E_hf_numerical(B, atom):
    # convert to Tesla
    B_tmp = B * 1e-4

    # make B iterable if it isnt already
    if not isinstance(B_tmp, Iterable):
        B_tmp = [B_tmp]

    I_id = qutip.identity(int(2*atom.I+1))
    J_id = qutip.identity(int(2*atom.J+1))

    # define operators for each spatial direction
    I_op = qutip.jmat(atom.I)

    I_x = qutip.tensor(I_op[0], qutip.identity(int(2*atom.J+1)))
    I_y = qutip.tensor(I_op[1], qutip.identity(int(2*atom.J+1)))
    I_z = qutip.tensor(I_op[2], qutip.identity(int(2*atom.J+1)))

    J_op = qutip.jmat(atom.J)

    J_x = qutip.tensor(qutip.identity(int(2*atom.I+1)), J_op[0])
    J_y = qutip.tensor(qutip.identity(int(2*atom.I+1)), J_op[1])
    J_z = qutip.tensor(qutip.identity(int(2*atom.I+1)), J_op[2])

    # F_x = J_x + I_x
    # F_y = J_y + I_y
    # F_z = J_z + I_z

    I2 = I_x**2 + I_y**2 + I_z**2
    J2 = J_x**2 + J_y**2 + J_z**2

    IJ = I_x*J_x+I_y*J_y+I_z*J_z


    # figure out the number of states
    F_s = np.arange(atom.I - atom.J, atom.I + atom.J + 1)
    n_states = sum([int(2*F+1) for F in F_s])

    E = np.zeros((n_states, len(B_tmp)))

    for i, b in enumerate(B_tmp):
        H_hf = h*atom.a_hf/hbar * (IJ)
        if atom.J != 1/2:
            H_hf += h*atom.b_hf/hbar * ( 3*(I2+J2+2*IJ) + 3/2*(IJ) - IJ**2) / (2*atom.I*(2*atom.I-1) * atom.J*(2*atom.J-1))

        H_z = (mu_bohr/hbar)*(qutip.tensor(I_id,atom.g_J*J_op[2]) + qutip.tensor(atom.g_I*I_op[2], J_id)) * b
        H = H_hf + H_z
        E[:, i] = H.eigenenergies()*hbar
    return E


def correct_traces(data, limit=None):
    n_traces = data.shape[0]
    diff2 = np.diff(data, n=2, axis=1)
    
    # set default limit (1 sigma)
    if not limit:
        limit = np.abs(diff2).mean() + 1*diff2.std()
    
    # identify outsiders
    ax_upper, index_upper = np.where(diff2>limit)
    ax_lower, index_lower = np.where(diff2<-limit)

    # make a dict that contains the traces to switch for each index
    switch_dict = {}
    for i, index in enumerate(index_upper):
        switch_dict[index] = [ax_upper[i]]
    for i, index in enumerate(index_lower):
        switch_dict[index].append(ax_lower[i])
    
    # keep track of the current trace identities
    trace_id_dict = {}
    for i in range(n_traces):
        trace_id_dict[i] = i
    
    new_data = data.copy()
    last_index = 0
    
    for index in sorted(switch_dict.keys()):
        index_sw = index + 1
        ax1, ax2 = switch_dict[index]
        if (index+1 in switch_dict.keys() and switch_dict[index+1] != switch_dict[index]) \
        or (index+1 not in switch_dict.keys()):
            tmp = new_data[trace_id_dict[ax1], index_sw:].copy()
            new_data[trace_id_dict[ax1], index_sw:] = new_data[trace_id_dict[ax2], index_sw:].copy()
            new_data[trace_id_dict[ax2], index_sw:] = tmp
#             trace_id_dict[trace_id_dict[ax1]] = trace_id_dict[ax2]
#             trace_id_dict[trace_id_dict[ax2]] = trace_id_dict[ax1]
            tmp_ax = trace_id_dict[ax1]
            trace_id_dict[ax1] = trace_id_dict[ax2]
            trace_id_dict[ax2] = tmp_ax
            
        last_index = index
    
    return new_data