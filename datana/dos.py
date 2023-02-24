# -*- coding: utf-8 -*-

import numpy as np

class Dos:
    def __init__(self):
        self.n_edos = 0
        self.e_fermi = 0.0
    
    def read_DOSCAR(self, file_name = "DOSCAR"):
        with open(file_name) as f:
            list_energy = []
            list_dos = []
            list_tdos = []
            for i in range(5):
                f.readline()
            line = f.readline().strip().split()
            self.n_edos = int(line[2])
            self.e_fermi = float(line[3])
            for i in range(self.n_edos):
                line = f.readline()
                list_energy.append(float(line.strip().split()[0]))
                list_dos.append(float(line.strip().split()[1]))
                list_tdos.append(float(line.strip().split()[2]))
            self._energy = np.array(list_energy) - self.e_fermi
            self._dos = np.array(list_dos)
            self._tdos = np.array(list_tdos)
            
    def get_dos(self):
        return np.stack((self._energy, self._dos), axis = 1)
    
    def get_tdos(self):
        return np.stack((self._energy, self._tdos), axis = 1)
