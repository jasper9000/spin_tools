e = 1.602176634e-19 # C
m_e = 9.1093837015e-31 #keg
m_p = 1.67262192369e-27 # kg
h = 6.62607015e-34 # m^2 kg / s
hbar = 1.0545718e-34 # m^2 kg / s
mu_bohr = e*hbar/(2*m_e)
mu_nuclear = e*hbar/(2*m_p)
g_L = 0.99998627
g_e = 2.0023193043737

def get_g_J(atom):
    return g_L * (atom.J*(atom.J+1) - atom.S*(atom.S+1) + \
                    atom.L*(atom.L+1))/(2*atom.J*(atom.J+1)) + \
                    g_e * (atom.J*(atom.J+1) + atom.S*(atom.S+1) - \
                    atom.L*(atom.L+1))/(2*atom.J*(atom.J+1))

class K40:
    def __init__(self):
        self.S = 1/2
        self.I = 4
        self.g_I = 0.000176490

class K40_4S_J12(K40):
    def __init__(self):
        super().__init__()
        self.L = 0
        self.J = 1/2
        self.g_J = get_g_J(self)

        self.a_hf = -285.7308*1e6 # h * Hz
        self.delta_E_hf = self.a_hf * (self.I+0.5) # Hz
        self.b_hf = 0

# class K40_4S_J12_F92(K40):
#     def __init__(self):
#         super().__init__()
#         self.L = 0
#         # self.S = 1/2
#         self.J = 1/2
#         # self.I = 4
#         self.F = 9/2

#         # self.g_I = 0.000176490
#         self.g_J = get_g_J(self)

#         self.a_hf = -285.7308*1e6 # h * Hz
#         self.delta_E_hf = self.a_hf * (self.I+0.5) # Hz
#         self.b_hf = 0

# class K40_4S_J12_F72(K40):
#     def __init__(self):
#         super().__init__()
#         self.L = 0
#         # self.S = 1/2
#         self.J = 1/2
#         # self.I = 4
#         self.F = 7/2

#         # self.g_I = 0.000176490
#         self.g_J = get_g_J(self)

#         self.a_hf = -285.7308*1e6 # h * Hz
#         self.delta_E_hf = self.a_hf * (self.I+0.5) # Hz
#         self.b_hf = 0


class K40_4P_J32(K40):
    def __init__(self):
        super().__init__()
        self.L = 1
        self.J = 3/2
        self.g_J = get_g_J(self)

        self.a_hf = -7.585*1e6
        self.b_hf = -3.445*1e6


########################################################################################

class Rb87:
    def __init__(self):
        self.S = 1/2
        self.I = 3/2
        self.g_I = -0.0009951414

class Rb87_5S_J12(Rb87):
    def __init__(self):
        super().__init__()
        self.L = 0
        self.J = 1/2
        self.g_J = get_g_J(self)

        self.a_hf = 3417.34130545215*1e6 # h * Hz
        self.delta_E_hf = self.a_hf * (self.I+0.5) # Hz
        self.b_hf = 0

# class Rb87_5S_J12_F1(Rb87):
#     def __init__(self):
#         super().__init__()
#         self.L = 0
#         self.J = 1/2
#         self.F = 1

#         # self.g_I = 0.000176490
#         self.g_J = g_L * (self.J*(self.J+1) - self.S*(self.S+1) + \
#                     self.L*(self.L+1))/(2*self.J*(self.J+1)) + \
#                     g_e * (self.J*(self.J+1) + self.S*(self.S+1) - \
#                     self.L*(self.L+1))/(2*self.J*(self.J+1))
#         self.g_F = -1/2

#         self.a_hf = 3417.34130545215*1e6 # h * Hz
#         self.delta_E_hf = self.a_hf * (self.I+0.5) # Hz
#         self.b_hf = 0

# class Rb87_5S_J12_F2(Rb87):
#     def __init__(self):
#         super().__init__()
#         self.L = 0
#         self.J = 1/2
#         self.F = 2

#         # self.g_I = 0.000176490
#         self.g_J = g_L * (self.J*(self.J+1) - self.S*(self.S+1) + \
#                     self.L*(self.L+1))/(2*self.J*(self.J+1)) + \
#                     g_e * (self.J*(self.J+1) + self.S*(self.S+1) - \
#                     self.L*(self.L+1))/(2*self.J*(self.J+1))
#         self.g_F = 1/2

#         self.a_hf = 3417.34130545215*1e6 # h
#         self.delta_E_hf = self.a_hf * (self.I+0.5) # Hz
#         self.b_hf = 0


class Rb87_5P_J32(Rb87):
    def __init__(self):
        super().__init__()
        self.L = 1
        self.J = 3/2
        self.g_J = get_g_J(self)
        
        self.a_hf = 84.7185*1e6 # h * Hz
        self.delta_E_hf = self.a_hf * (self.I+0.5) * 1e6 # Hz
        self.b_hf = 12.4965*1e6 # h * Hz