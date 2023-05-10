# -*- coding: utf-8 -*-

import numpy as np
import copy
from datana.structure import Structure

class Dynamic:
    def __init__(self):
        self._temperature = []
        self._energy = []
        self._structures = []
        
    def read_OSZICAR_XDATCAR(self, OSZICAR = "OSZICAR", XDATCAR = "XDATCAR"):
        list_temperature = []
        list_energy = []
        list_structure = []
        if(OSZICAR):
            with open(OSZICAR) as f:
                for line in f:
                    if("T=" in line):
                        list_line = line.strip().split()
                        list_temperature.append(float(list_line[2]))
                        list_energy.append(float(list_line[4]))
        else:
            pass
        if(XDATCAR):
            with open(XDATCAR) as f:
                s_single = Structure()
                temp_lattice = []
                for i in range(2):
                    f.readline()
                for i in range(3):
                    data = f.readline().strip().split()
                    data = list(map(float, data))
                    temp_lattice.append(data)
                s_single._lattice = np.array(temp_lattice)
                list_species = f.readline().strip().split()
                num_atom = f.readline().strip().split()
                num_atom = list(map(int, num_atom))
                for index, value in enumerate(list_species):
                    s_single._species += [value] * num_atom[index]
                for line in f:
                    temp_coords_direct = []
                    for i in range(int(np.sum(num_atom))):
                        data = f.readline().split()
                        temp_coords_direct.append(data)
                    s_single._coords_direct = np.array(temp_coords_direct).astype('float64')
                    s_single.update_cartesian()
                    s_single_copy = copy.deepcopy(s_single)
                    list_structure.append(s_single_copy)
        else:
            pass
        if(not len(list_temperature)):
            list_temperature, list_energy = [0] * len(list_structure), [0] * len(list_structure)
        elif(not len(list_structure)):
            s_single = Structure()
            list_structure = [s_single] * len(list_temperature)
        if(len(list_energy) != len(list_structure)):
            print("ERROR: Energy and structure list must the same size...")
            return        
        self._temperature = np.array(list_temperature)
        self._energy = np.array(list_energy)
        self._structures = list_structure
    
    def update_OSZICAR(self, OSZICAR = "OSZICAR"):
        list_temperature = []
        list_energy = []
        with open(OSZICAR) as f:
            for line in f:
                if("T=" in line):
                    list_line = line.strip().split()
                    list_temperature.append(float(list_line[2]))
                    list_energy.append(float(list_line[4]))
        if(len(list_energy) == len(self._structures)):
            self._temperature = np.array(list_temperature)
            self._energy = np.array(list_energy)
    
    def update_XDATCAR(self, XDATCAR = "XDATCAR"):
        list_structure = []
        with open(XDATCAR) as f:
            s_single = Structure()
            temp_lattice = []
            for i in range(2):
                f.readline()
            for i in range(3):
                data = f.readline().strip().split()
                data = list(map(float, data))
                temp_lattice.append(data)
            s_single._lattice = np.array(temp_lattice)
            list_species = f.readline().strip().split()
            num_atom = f.readline().strip().split()
            num_atom = list(map(int, num_atom))
            for index, value in enumerate(list_species):
                s_single._species += [value] * num_atom[index]
            for line in f:
                temp_coords_direct = []
                for i in range(int(np.sum(num_atom))):
                    data = f.readline().split()
                    temp_coords_direct.append(data)
                s_single._coords_direct = np.array(temp_coords_direct).astype('float64')
                s_single.update_cartesian()
                s_single_copy = copy.deepcopy(s_single)
                list_structure.append(s_single_copy)
        if(len(self._structures) == len(list_energy)):
            self._structures = list_structure
    
    def get_temperature(self):
        return self._temperature
        
    def get_energy(self):
        return self._energy
    
    def cal_msd(self, index_atom = -1):
        num_atom = len(self._structures[0]._species)
        res_msd = []
        for item in self._structures:
            if(len(item._species) != num_atom):
                print("ERROR: cannot calculate msd with different atom number...")
                break
            if(index_atom >= 0 and index_atom < num_atom):
                x2 = (item._coords_cartesian[index_atom][0] - self._structures[0]._coords_cartesian[index_atom][0])**2
                y2 = (item._coords_cartesian[index_atom][1] - self._structures[0]._coords_cartesian[index_atom][1])**2
                z2 = (item._coords_cartesian[index_atom][2] - self._structures[0]._coords_cartesian[index_atom][2])**2
                res_msd.append(np.sqrt(x2+y2+z2))
            else:
                msd_single = []
                for index, coords_atom in enumerate(item._coords_cartesian):
                    x2 = (coords_atom[0] - self._structures[0]._coords_cartesian[index][0]) ** 2
                    y2 = (coords_atom[1] - self._structures[0]._coords_cartesian[index][1]) ** 2
                    z2 = (coords_atom[2] - self._structures[0]._coords_cartesian[index][2]) ** 2
                    msd_single.append(np.sqrt(x2 + y2 + z2))
                res_msd.append(np.mean(msd_single))
        return np.array(res_msd)