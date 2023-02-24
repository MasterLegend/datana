# -*- coding: utf-8 -*-

# 原子半径参考： E. Clementi, D.L.Raimondi, 和 W.P. Reinhardt, 《化学物理期刊》（J. Chem. Phys） 1963, 38, 2686.
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
    "Cs":2.98, "Ba":2.53
}

class Atom:
    def __init__(self, atom_name = None):
        self._type = atom_name
        self._radius = map_atom_radius.get(atom_name, 0.0)
    
    def __str__(self):
        str_summary = 'Atom Summary\n'
        str_summary += f'Type:\t{self._type}\n'
        str_summary += f"Radius:\t{self._radius}\n"
        return str_summary
    
    def __repr__(self):
        return self.__str__()
    
    def get_type(self):
        return self._type
    
    def get_radius(self):
        return self._radius