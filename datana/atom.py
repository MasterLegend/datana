# -*- coding: utf-8 -*-

import numpy as np

map_atom_number = {
    "H":1, "He":2,
    "Li":3, "Be":4, "B":5, "C":6, "N":7, "O":8, "F":9, "Ne":10,
    "Na":11, "Mg":12, "Al":13, "Si":14, "P":15, "S":16, "Cl":17, "Ar":18,
    "K":19, "Ca":20,
    "Sc":21, "Ti":22, "V":23, "Cr":24, "Mn":25, "Fe":26, "Co":27, "Ni":28, "Cu":29, "Zn":30,
    "Ga":31, "Ge":32, "As":33, "Se":34, "Br":35, "Kr":36,
    "Rb":37, "Sr":38,
    "Y":39, "Zr":40, "Nb":41, "Mo":42, "Tc":43, "Ru":44, "Rh":45, "Pd":46, "Ag":47, "Cd":48,
    "In":49, "Sn":50, "Sb":51, "Te":52, "I":53, "Xe":54,
    "Cs":55, "Ba":56,
    "La":57, "Ce":58, "Pr":59, "Nd":60, "Pm":61, "Sm":62, "Eu":63, "Gd":64,
    "Tb":65, "Dy":66, "Ho":67, "Er":68, "Tm":69, "Yb":70, "Lu":71,
    "Hf":72, "Ta":73, "W":74, "Re":75, "Os":76, "Ir":77, "Pt":78, "Au":79, "Hg":80,
    "Tl":81, "Pb":82, "Bi":83,
}

map_atom_periodic_number = {
    "H":1, "He":1,
    "Li":2, "Be":2, "B":2, "C":2, "N":2, "O":2, "F":2, "Ne":2,
    "Na":3, "Mg":3, "Al":3, "Si":3, "P":3, "S":3, "Cl":3, "Ar":3,
    "K":4, "Ca":4,
    "Sc":4, "Ti":4, "V":4, "Cr":4, "Mn":4, "Fe":4, "Co":4, "Ni":4, "Cu":4, "Zn":4,
    "Ga":4, "Ge":4, "As":4, "Se":4, "Br":4, "Kr":4,
    "Rb":5, "Sr":5,
    "Y":5, "Zr":5, "Nb":5, "Mo":5, "Tc":5, "Ru":5, "Rh":5, "Pd":5, "Ag":5, "Cd":5,
    "In":5, "Sn":5, "Sb":5, "Te":5, "I":5, "Xe":5,
    "Cs":6, "Ba":6,
    "La":6, "Ce":6, "Pr":6, "Nd":6, "Pm":6, "Sm":6, "Eu":6, "Gd":6,
    "Tb":6, "Dy":6, "Ho":6, "Er":6, "Tm":6, "Yb":6, "Lu":6,
    "Hf":6, "Ta":6, "W":6, "Re":6, "Os":6, "Ir":6, "Pt":6, "Au":6, "Hg":6,
    "Tl":6, "Pb":6, "Bi":6,
}

map_atom_group_number = {
    "H":1, "He":18,
    "Li":1, "Be":2, "B":13, "C":14, "N":15, "O":16, "F":17, "Ne":18,
    "Na":1, "Mg":2, "Al":13, "Si":14, "P":15, "S":16, "Cl":17, "Ar":18,
    "K":1, "Ca":2,
    "Sc":3, "Ti":4, "V":5, "Cr":6, "Mn":7, "Fe":8, "Co":9, "Ni":10, "Cu":11, "Zn":12,
    "Ga":13, "Ge":14, "As":15, "Se":16, "Br":17, "Kr":18,
    "Rb":1, "Sr":2,
    "Y":3, "Zr":4, "Nb":5, "Mo":6, "Tc":7, "Ru":8, "Rh":9, "Pd":10, "Ag":11, "Cd":12,
    "In":13, "Sn":14, "Sb":15, "Te":16, "I":17, "Xe":18,
    "Cs":1, "Ba":2,
    "La":3, "Ce":3, "Pr":3, "Nd":3, "Pm":3, "Sm":3, "Eu":3, "Gd":3,
    "Tb":3, "Dy":3, "Ho":3, "Er":3, "Tm":3, "Yb":3, "Lu":3,
    "Hf":4, "Ta":5, "W":6, "Re":7, "Os":8, "Ir":9, "Pt":10, "Au":11, "Hg":12,
    "Tl":13, "Pb":14, "Bi":15,
}

map_atom_weight = {
    "H":1.00794, "He":4.002602,
    "Li":6.941, "Be":9.012182, "B":10.811, "C":12.0107, "N":14.0067, "O":15.9994, "F":18.9984032, "Ne":20.1797,
    "Na":22.989770, "Mg":24.3050, "Al":26.981538, "Si":28.0855, "P":30.973761, "S":32.065, "Cl":35.453, "Ar":39.948,
    "K":39.0983, "Ca":40.078,
    "Sc":44.955910, "Ti":47.867, "V":50.9415, "Cr":51.9961, "Mn":54.938049, "Fe":55.845, "Co":58.933200, "Ni":58.6934, "Cu":63.546, "Zn":65.409,
    "Ga":69.723, "Ge":72.64, "As":74.92160, "Se":78.96, "Br":79.904, "Kr":83.798,
    "Rb":85.4678, "Sr":87.62,
    "Y":88.90585, "Zr":91.224, "Nb":92.90638, "Mo":95.94, "Tc":98, "Ru":101.07, "Rh":102.90550, "Pd":106.42, "Ag":107.8682, "Cd":112.411,
    "In":114.818, "Sn":118.710, "Sb":121.760, "Te":127.60, "I":126.90447, "Xe":131.293,
    "Cs":132.90545, "Ba":137.327,
    "La":138.9055, "Ce":140.116, "Pr":140.90765, "Nd":144.24, "Pm":145, "Sm":150.36, "Eu":151.964, "Gd":157.25,
    "Tb":158.92534, "Dy":162.500, "Ho":164.93032, "Er":167.259, "Tm":168.93421, "Yb":173.04, "Lu":174.967,
    "Hf":178.49, "Ta":180.9479, "W":183.84, "Re":186.207, "Os":190.23, "Ir":192.217, "Pt":195.078, "Au":196.96655, "Hg":200.59,
    "Tl":204.3833, "Pb":207.2, "Bi":208.98038,
}

map_atom_radius = {
    "H":0.53, "He":0.31,
    "Li":1.67, "Be":1.12, "B":0.87, "C":0.67, "N":0.56, "O":0.48, "F":0.42, "Ne":0.38,
    "Na":1.90, "Mg":1.45, "Al":1.18, "Si":1.11, "P":0.98, "S":0.88, "Cl":0.79, "Ar":0.71,
    "K":2.43, "Ca":1.94,
    "Sc":1.84, "Ti":1.76, "V":1.71, "Cr":1.66, "Mn":1.61, "Fe":1.56, "Co":1.52, "Ni":1.49, "Cu":1.45, "Zn":1.42,
    "Ga":1.36, "Ge":1.25, "As":1.14, "Se":1.03, "Br":0.94, "Kr":0.88,
    "Rb":2.65, "Sr":2.19,
    "Y":2.12, "Zr":2.06, "Nb":1.98, "Mo":1.90, "Tc":1.83, "Ru":1.78, "Rh":1.73, "Pd":1.69, "Ag":1.65, "Cd":1.61,
    "In":1.56, "Sn":1.45, "Sb":1.33, "Te":1.23, "I":1.15, "Xe":1.08,
    "Cs":2.98, "Ba":2.53,
    "La":2.26, "Ce":2.10, "Pr":2.47, "Nd":2.06, "Pm":2.05, "Sm":2.38, "Eu":2.31, "Gd":2.33,
    "Tb":2.25, "Dy":2.28, "Ho":2.26, "Er":2.26, "Tm":2.22, "Yb":2.22, "Lu":2.17,
    "Hf":2.08, "Ta":2.00, "W":1.93, "Re":1.88, "Os":1.85, "Ir":1.80, "Pt":1.77, "Au":1.74, "Hg":1.71,
    "Tl":1.56, "Pb":1.54, "Bi":1.43,
}

map_atom_electron_affinity = {
    "H":0.754195, "He":np.nan,
    "Li":0.618049, "Be":np.nan, "B":0.279723, "C":1.262119, "N":np.nan, "O":1.4611096, "F":3.4011895, "Ne":np.nan,
    "Na":0.547926, "Mg":np.nan, "Al":0.43283, "Si":1.3895220, "P":0.7465, "S":2.077103, "Cl":3.612724, "Ar":np.nan,
    "K":0.50147, "Ca":0.02455,
    "Sc":0.188, "Ti":0.079, "V":0.525, "Cr":0.666, "Mn":np.nan, "Fe":0.151, "Co":0.662, "Ni":1.156, "Cu":1.235, "Zn":np.nan,
    "Ga":0.43, "Ge":1.232712, "As":0.814, "Se":2.020670, "Br":3.363588, "Kr":np.nan,
    "Rb":0.48592, "Sr":0.048,
    "Y":0.307, "Zr":0.426, "Nb":0.893, "Mo":0.748, "Tc":0.55, "Ru":1.05, "Rh":1.137, "Pd":0.562, "Ag":1.302, "Cd":np.nan,
    "In":0.3, "Sn":1.112067, "Sb":1.046, "Te":1.9708, "I":3.059037, "Xe":np.nan,
    "Cs":0.471626, "Ba":0.14462,
    "La":0.47, "Ce":np.nan, "Pr":np.nan, "Nd":np.nan, "Pm":np.nan, "Sm":np.nan, "Eu":np.nan, "Gd":np.nan,
    "Tb":np.nan, "Dy":np.nan, "Ho":np.nan, "Er":np.nan, "Tm":np.nan, "Yb":-0.020, "Lu":0.34,
    "Hf":0, "Ta":0.322, "W":0.815, "Re":0.15, "Os":1.1, "Ir":1.5638, "Pt":2.128, "Au":2.30863, "Hg":np.nan,
    "Tl":0.2, "Pb":0.364, "Bi":0.946,
}

map_atom_ionization_potential = {
    "H":13.59844, "He":24.58741,
    "Li":5.39172, "Be":9.3227, "B":8.29803, "C":11.26030, "N":14.53414, "O":13.61806, "F":17.42282, "Ne":21.5646,
    "Na":5.13908, "Mg":7.64624, "Al":5.98577, "Si":8.15169, "P":10.48669, "S":10.36001, "Cl":12.96764, "Ar":15.75962,
    "K":4.34066, "Ca":6.11316,
    "Sc":6.5615, "Ti":6.8281, "V":6.7462, "Cr":6.7665, "Mn":7.43402, "Fe":7.9024, "Co":7.8810, "Ni":7.6398, "Cu":7.72638, "Zn":9.3942,
    "Ga":5.99930, "Ge":7.8994, "As":9.7886, "Se":9.75238, "Br":11.81381, "Kr":13.99961,
    "Rb":4.17713, "Sr":5.6949,
    "Y":6.2171, "Zr":6.63390, "Nb":6.75885, "Mo":7.09243, "Tc":7.28, "Ru":7.36050, "Rh":7.45890, "Pd":8.3369, "Ag":7.5762, "Cd":8.9938,
    "In":5.78636, "Sn":7.3439, "Sb":8.6084, "Te":9.0096, "I":10.45126, "Xe":12.1298,
    "Cs":3.89390, "Ba":5.21170,
    "La":5.5769, "Ce":5.5387, "Pr":5.473, "Nd":5.5250, "Pm":5.582, "Sm":5.6436, "Eu":5.6704, "Gd":6.1501,
    "Tb":5.8638, "Dy":5.9389, "Ho":6.0215, "Er":6.1077, "Tm":6.18431, "Yb":6.25416, "Lu":5.4259,
    "Hf":6.82507, "Ta":7.5496, "W":7.8640, "Re":7.8335, "Os":8.4382, "Ir":8.9670, "Pt":8.9587, "Au":9.2255, "Hg":10.43750,
    "Tl":6.1082, "Pb":7.41666, "Bi":7.2856,
}

map_atom_polarizability = {
    "H":0.666793, "He":0.204956,
    "Li":24.3, "Be":5.60, "B":3.03, "C":1.76, "N":1.10, "O":0.802, "F":0.557, "Ne":0.3956,
    "Na":24.11, "Mg":10.6, "Al":6.8, "Si":5.38, "P":3.63, "S":2.90, "Cl":2.18, "Ar":1.6411,
    "K":43.4, "Ca":22.8,
    "Sc":17.8, "Ti":14.6, "V":12.4, "Cr":11.6, "Mn":9.4, "Fe":8.4, "Co":7.5, "Ni":6.8, "Cu":6.2, "Zn":5.75,
    "Ga":8.12, "Ge":6.07, "As":4.31, "Se":3.77, "Br":3.05, "Kr":2.4844,
    "Rb":47.3, "Sr":27.6,
    "Y":22.7, "Zr":17.9, "Nb":15.7, "Mo":12.8, "Tc":11.4, "Ru":9.6, "Rh":8.6, "Pd":4.8, "Ag":7.2, "Cd":7.36,
    "In":10.2, "Sn":7.7, "Sb":6.6, "Te":5.5, "I":5.35, "Xe":4.044,
    "Cs":59.42, "Ba":39.7,
    "La":31.1, "Ce":29.6, "Pr":28.2, "Nd":31.4, "Pm":30.1, "Sm":28.8, "Eu":27.7, "Gd":23.5,
    "Tb":25.5, "Dy":24.5, "Ho":23.6, "Er":22.7, "Tm":21.8, "Yb":21.0, "Lu":21.9,
    "Hf":16.2, "Ta":13.1, "W":11.1, "Re":9.7, "Os":8.5, "Ir":7.6, "Pt":6.5, "Au":5.8, "Hg":5.02,
    "Tl":7.6, "Pb":6.8, "Bi":7.4,
}

class Atom:
    def __init__(self, atom_name = None):
        self.__type = atom_name
        self.__number = map_atom_number.get(atom_name, 0)
        self.__periodic = map_atom_periodic_number.get(atom_name, 0)
        self.__group = map_atom_group_number.get(atom_name, 0)
        self.__weight = map_atom_weight.get(atom_name, 0.0)
        self.__radius = map_atom_radius.get(atom_name, 0.0)
        self.__electron_affinity = map_atom_electron_affinity.get(atom_name, np.nan)
        self.__ionization_potential = map_atom_ionization_potential.get(atom_name, 0.0)
        self.__polarizability = map_atom_polarizability.get(atom_name, 0.0)
    
    def __str__(self):
        str_summary = 'Atom Summary\n'
        str_summary += f'Type:\t{self.__number}-{self.__type}\n'
        str_summary += f"Periodic: {self.__periodic}\tGroup: {self.__group}\n"
        str_summary += f"Weight: {self.__weight}\tRadius: {self.__radius}\n"
        str_summary += f"EA: {self.__electron_affinity}\tIP: {self.__ionization_potential}\tPolarizability: {self.__polarizability}\n"
        return str_summary
    
    def __repr__(self):
        return self.__str__()
    
    def get_type(self):
        return self.__type
    
    def get_number(self):
        return self.__number
    
    def get_periodic(self):
        return self.__periodic
    
    def get_group(self):
        return self.__group
    
    def get_weight(self):
        return self.__weight
    
    def get_radius(self):
        return self.__radius
    
    def get_electron_affinity(self):
        return self.__electron_affinity
    
    def get_ionization_potential(self):
        return self.__ionization_potential
    
    def get_polarizability(self):
        return self.__polarizability
