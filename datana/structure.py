# -*- coding: utf-8 -*-

import numpy as np
from functools import reduce
import copy
from datana.atom import Atom
import random

class Structure:
    def __init__(self, a = 1.0, b = 1.0, c = 1.0, alpha = np.pi/2, beta = np.pi/2, gamma = np.pi/2):
        self._lattice = np.array([[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]])
        self._atom = []
        self._species = []
        self._coords_direct = []
        self._coords_cartesian = []
        self._forces = []
        self._energy = 0.0
        self._isFix = []
    
    def __str__(self):
        str_summary = 'Structure Summary\n'
        str_summary += 'Lattice\n'
        abc = self.get_abc()
        str_summary += '  abc:\t'
        str_summary += f'{abc[0]:.10f}\t{abc[1]:.10f}\t{abc[2]:.10f}\n'
        angles = self.get_angles()
        str_summary += '  ang:\t'
        str_summary += f'{angles[0]:.10f}\t{angles[1]:.10f}\t{angles[2]:.10f}\n'
        str_summary += '  vec_a:\t'
        str_summary += f'{self._lattice[0][0]:.10f}\t{self._lattice[0][1]:.10f}\t{self._lattice[0][2]:.10f}\n'
        str_summary += '  vec_b:\t'
        str_summary += f'{self._lattice[1][0]:.10f}\t{self._lattice[1][1]:.10f}\t{self._lattice[1][2]:.10f}\n'
        str_summary += '  vec_c:\t'
        str_summary += f'{self._lattice[2][0]:.10f}\t{self._lattice[2][1]:.10f}\t{self._lattice[2][2]:.10f}\n'
        str_summary += f'Atom  {len(self._species)} atom(s) in total\n'
        for i in range(len(self._species)):
            str_summary += f'  {self._species[i]}\t'
            str_summary += f'{self._coords_direct[i][0]:.10f}\t{self._coords_direct[i][1]:.10f}\t{self._coords_direct[i][2]:.10f}\n'
        return str_summary
    
    def __repr__(self):
        return self.__str__()
    
    def get_atom(self, i):
        if(i < len(self._atom)):
            return self._atom[i]
    
    def get_abc(self):
        return np.sqrt((self._lattice * self._lattice).sum(axis = 1))
    
    def get_angles(self):
        abc = np.sqrt((self._lattice * self._lattice).sum(axis = 1))
        alpha = np.arccos(self._lattice[1].dot(self._lattice[2])/(abc[1] * abc[2]))
        beta = np.arccos(self._lattice[0].dot(self._lattice[2])/(abc[0] * abc[2]))
        gamma = np.arccos(self._lattice[0].dot(self._lattice[1])/(abc[0] * abc[1]))
        return np.array([alpha, beta, gamma]) * 180 / np.pi
    
    def read_POSCAR(self, file_name = 'POSCAR'):
        temp_lattice = []
        temp_coords_direct = []
        temp_coords_cartesian = []
        isFix = False
        
        with open(file_name) as f:
            for i in range(2):
                f.readline()
            for i in range(3):
                data = f.readline().strip().split()
                data = list(map(float, data))
                temp_lattice.append(data)
            self._lattice = np.array(temp_lattice)
            list_species = f.readline().strip().split()
            num_atom = f.readline().strip().split()
            num_atom = list(map(int, num_atom))
            for index, value in enumerate(list_species):
                self._species += [value] * num_atom[index]
                self._atom += [Atom(value)] * num_atom[index]  
            isDirect = f.readline().strip()
            if(isDirect[0] == 'S'):
                isFix = True
                isDirect = f.readline().strip()
            for i in range(int(np.sum(num_atom))):
                data = f.readline().split()
                if(isFix):
                    if(isDirect[0] in ['C', 'c', 'K', 'k']):
                        temp_coords_cartesian.append(data[:3])
                    else:
                        temp_coords_direct.append(data[:3])
                    self._isFix.append(data[3:])
                else:
                    if(isDirect[0] in ['C', 'c', 'K', 'k']):
                        temp_coords_cartesian.append(data)
                    else:
                        temp_coords_direct.append(data)
                    self._isFix.append(['T','T','T'])
            if(isDirect[0] in ['C', 'c', 'K', 'k']):
                self._coords_cartesian = np.array(temp_coords_cartesian).astype('float64')
            else:
                self._coords_direct = np.array(temp_coords_direct).astype('float64')
        self.fillna()
    
    def fillna(self):
        if(len(self._coords_cartesian) > len(self._coords_direct)):
            self._coords_direct = self._coords_cartesian @ np.linalg.inv(self._lattice)
        elif(len(self._coords_direct) > len(self._coords_cartesian)):
            self._coords_cartesian = self._coords_direct @ self._lattice
        if(len(self._forces) != len(self._species)):
            self._forces = np.zeros(self._coords_direct.shape)
        if(len(self._isFix) != len(self._species)):
            self._isFix = np.array([['T','T','T'] for i in range(len(self._coords_direct))])
    
    def update_cartesian(self):
        self._coords_cartesian = self._coords_direct @ self._lattice
    
    def update_direct(self):
        self._coords_direct = self._coords_cartesian @ np.linalg.inv(self._lattice)
    
    def cal_dist(self, i = 0, j = 0):
        num_atom = len(self._coords_cartesian)
        if(i < 0 or j < 0 or i >= num_atom or j >= num_atom):
            print('Error: atom number out of index...')
            return
        else:
            x2 = (self._coords_cartesian[i][0] - self._coords_cartesian[j][0]) * (self._coords_cartesian[i][0] - self._coords_cartesian[j][0])
            y2 = (self._coords_cartesian[i][1] - self._coords_cartesian[j][1]) * (self._coords_cartesian[i][1] - self._coords_cartesian[j][1])
            z2 = (self._coords_cartesian[i][2] - self._coords_cartesian[j][2]) * (self._coords_cartesian[i][2] - self._coords_cartesian[j][2])
            return np.sqrt(x2 + y2 + z2)
    
    def cal_volume(self):
        vol = np.sum(np.cross(self._lattice[0], self._lattice[1]) * self._lattice[2])
        return abs(vol)
    
    def re_orientation(self):
        pass
    
    def incell(self):
        while((self._coords_direct > 1).any()):
            self._coords_direct[np.where(self._coords_direct > 1.0)] -= 1
        while((self._coords_direct < 0).any()):
            self._coords_direct[np.where(self._coords_direct < 1.0)] += 1
        self.update_cartesian()
    
    def fix(self, threshold = 0.0, direction = 'z'):
        for i in range(len(self._coords_direct)):
            if(self._coords_direct[i][2] < threshold):
                self._isFix[i] = ['F','F','F']
            else:
                self._isFix[i] = ['T','T','T']
    
    def supercell(self, a = 1, b = 1, c = 1):  # do not finish
        if(a <= 0 or b <= 0 or c <= 0):
            print("ERROR: positive multiples are required")
            return
        a, b, c = int(a), int(b), int(c)
        self._lattice[0] *= a
        self._lattice[1] *= b
        self._lattice[2] *= c
        self._species *= a*b*c
        self._coords_direct[:, 0] /= a
        self._coords_direct[:, 1] /= b
        self._coords_direct[:, 2] /= c
        coords_copy = copy.deepcopy(self._coords_direct)
        for i in range(a):
            for j in range(b):
                for k in range(c):
                    if(i or j or k):
                        self._coords_direct = np.concatenate((self._coords_direct, coords_copy + (i/a, j/b, k/c)))
        del coords_copy
        self.update_cartesian()
        force_copy = copy.deepcopy(self._forces)
        for i in range(a*b*c):
            self._forces = np.concatenate((self._forces, force_copy))
        del force_copy
        self._isFix *= a*b*c

    def replacement(self, initial_num = -1, replace_element = 'H'):
        self._species[initial_num] = replace_element
        self._atom[initial_num] = Atom(replace_element)
    
    def random_replacement(self, initial_element, replace_element, replace_num = -1):
        initial_element_num = 0
        replace_ele_list = []
        for i, j in enumerate(self._atom):
            if j._type == initial_element:
                initial_element_num += 1
                replace_ele_list.append(i)
        if replace_num < 0 or replace_num > initial_element_num: replace_number = initial_element_num
        replace_list = random.sample(replace_ele_list, replace_num)
        structures_from_random_replacement = Structures()
        for i in replace_list:
            structure_temp = copy.deepcopy(self)
            structure_temp.replacement(i, replace_element)
            structures_from_random_replacement.add(copy.deepcopy(structure_temp))
        return structures_from_random_replacement
    
    def to_POSCAR(self, file_name = 'POSCAR', mode = 'Direct'):
        isFix = False
        
        with open(file_name, 'w') as f:
            f.write('POSCAR_created_by_datana\n')
            f.write('1.0000000000\n')
            f.write(f'\t{self._lattice[0][0]:.10f}\t{self._lattice[0][1]:.10f}\t{self._lattice[0][2]:.10f}\n')
            f.write(f'\t{self._lattice[1][0]:.10f}\t{self._lattice[1][1]:.10f}\t{self._lattice[1][2]:.10f}\n')
            f.write(f'\t{self._lattice[2][0]:.10f}\t{self._lattice[2][1]:.10f}\t{self._lattice[2][2]:.10f}\n')
            reduce_species = reduce(lambda x, y: x if y in x else x + [y], [[], ] + self._species)
            for i in range(len(reduce_species)):
                f.write(f'  {reduce_species[i]}')
            f.write('\n')
            for i in range(len(reduce_species)):
                f.write(f'  {self._species.count(reduce_species[i])}')
            f.write('\n')
            if(np.sum(np.array(self._isFix) == 'F')):
                f.write('Selective dynamics\n')
                isFix = True
            f.write(f'{mode}\n')
            if(isFix):
                if(mode == 'Cartesian'):
                    for i in range(len(self._coords_cartesian)):
                        f.write(f'  {self._coords_cartesian[i][0]:.10f}  {self._coords_cartesian[i][1]:.10f}  {self._coords_cartesian[i][2]:.10f}  {self._isFix[i][0]}  {self._isFix[i][1]}  {self._isFix[i][2]}\n')
                else:
                    for i in range(len(self._coords_direct)):
                        f.write(f'  {self._coords_direct[i][0]:.10f}  {self._coords_direct[i][1]:.10f}  {self._coords_direct[i][2]:.10f}  {self._isFix[i][0]}  {self._isFix[i][1]}  {self._isFix[i][2]}\n')
            else:
                if(mode == 'Cartesian'):
                    for i in range(len(self._coords_cartesian)):
                        f.write(f'  {self._coords_cartesian[i][0]:.10f}  {self._coords_cartesian[i][1]:.10f}  {self._coords_cartesian[i][2]:.10f}\n')
                else:
                    for i in range(len(self._coords_direct)):
                        f.write(f'  {self._coords_direct[i][0]:.10f}  {self._coords_direct[i][1]:.10f}  {self._coords_direct[i][2]:.10f}\n')
    
    def to_cif(self, file_name = 'structure.cif'):
        with open(file_name, 'w') as f:
            f.write('data_created_by_datana\n')
            f.write('loop_\n')
            f.write('_symmetry_equiv_pos_as_xyz\n')
            f.write('  x,y,z\n')
            abc = self.get_abc()
            f.write(f'_cell_length_a\t{abc[0]:.4f}\n')
            f.write(f'_cell_length_b\t{abc[1]:.4f}\n')
            f.write(f'_cell_length_c\t{abc[2]:.4f}\n')
            angles = self.get_angles()
            f.write(f'_cell_angle_alpha\t{angles[0]:.4f}\n')
            f.write(f'_cell_angle_beta\t{angles[1]:.4f}\n')
            f.write(f'_cell_angle_gamma\t{angles[2]:.4f}\n')
            f.write('loop_\n')
            f.write('_atom_site_label\n')
            f.write('_atom_site_type_symbol\n')
            f.write('_atom_site_fract_x\n')
            f.write('_atom_site_fract_y\n')
            f.write('_atom_site_fract_z\n')
            f.write('_atom_site_U_iso_or_equiv\n')
            f.write('_atom_site_adp_type\n')
            f.write('_atom_site_occupancy\n')
            for i in range(len(self._coords_cartesian)):
                f.write(f'Atom{i}\t{self._species[i]}\t{self._coords_direct[i][0]:.5f}\t{self._coords_direct[i][1]:.5f}\t{self._coords_direct[i][2]:.5f}\t0.00000\tUiso\t1.00\n')

class Structures:
    def __init__(self):
        self._structures = []

    def __str__(self):
        str_summary = 'Structures Summary\n'
        str_summary += f'Frames:\t{len(self._structures)}'
        return str_summary
        
    def __repr__(self):
        return self.__str__()
    
    def add(self, struct):
        if(type(struct) == Structure):
            self._structures.append(struct)
        else:
            print('type input error...')
    
    def pop(self):
        try:
            self._structures.pop()
        except:
            print('pop error...')
    
    def read_OUTCAR(self, file_name = "OUTCAR", verbose = 1):
        with open(file_name) as f:
            s_single = Structure()
            list_species = []
            num_atom = [0]
            is_read_lattice_vector = False  # do not read in first time
            
            if(verbose):
                print('Reading the OUTCAR file...')
            for line in f:
                if("VRHFIN" in line):
                    list_species.append(line.split()[1].replace("=", "").replace(":", ""))
                elif("ions per type" in line):
                    num_atom = [int(i) for i in line.split()[4:]]
                if(sum(num_atom) > 0 and not s_single._species):
                    for index, value in enumerate(list_species):
                        s_single._species += [value] * num_atom[index]
                        s_single._atom += [Atom(value)] * num_atom[index]
                
                if("direct lattice vectors" in line):
                    if(not is_read_lattice_vector):
                        is_read_lattice_vector = True
                        continue
                    list_lattice = []
                    for i in range(3):
                        line = f.readline()
                        list_line = line.split()
                        list_lattice.append([float(list_line[0]), float(list_line[1]), float(list_line[2])])
                    s_single._lattice = np.array(list_lattice)
                elif("POSITION" in line and "TOTAL-FORCE" in line):
                    positions_single = []
                    forces_single = []
                    line = f.readline()
                    line = f.readline()
                    for i in range(sum(num_atom)):
                        list_line = line.strip('\n').split()
                        positions_single.append([float(list_line[0]), float(list_line[1]), float(list_line[2])])
                        forces_single.append([float(list_line[3]), float(list_line[4]), float(list_line[5])])
                        line = f.readline()
                    s_single._coords_cartesian = np.array(positions_single)
                    s_single._forces = np.array(forces_single)
                    for i in range(12):
                        line = f.readline()
                    s_single._energy = float(line.strip('\n').split()[-1])
                    s_single.fillna()
                    s_single.update_direct()
                    s_single_copy = copy.deepcopy(s_single)
                    self._structures.append(s_single_copy)
                    if(verbose):
                        print(f'\r{len(self._structures)}', end = '')
        if(verbose):
            print()
    
    def read_XDATCAR(self, file_name = "XDATCAR"):
        with open(file_name) as f:
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
                s_single._atom += [Atom(value)] * num_atom[index]
            for line in f:
                temp_coords_direct = []
                for i in range(int(np.sum(num_atom))):
                    data = f.readline().split()
                    temp_coords_direct.append(data)
                s_single._coords_direct = np.array(temp_coords_direct).astype('float64')
                s_single.update_cartesian()
                s_single_copy = copy.deepcopy(s_single)
                self._structures.append(s_single_copy)
    
    def read_xyz(self, file_name, verbose = 1):
        with open(file_name) as f:
            if(verbose):
                print('Reading the xyz file...')
            
            s_single = Structure()
            for line in f:
                num_atom = int(line.strip('\n'))
                line = f.readline()
                string_lattice = np.array(line.strip('\n').split()[1:])
                string_lattice[-1] = string_lattice[-1][:-1]     # remove the double quotation mark
                string_lattice = string_lattice.astype('float64')
                s_single._lattice = string_lattice.reshape(3, 3)
                positions_single = []
                for i in range(num_atom):
                    line = f.readline()
                    add_atom = line.strip('\n').split()
                    s_single._species.append(add_atom[0])
                    s_single._atom.append([Atom(add_atom[0])])
                    positions_single.append(list(map(float, add_atom[1:])))
                s_single._coords_cartesian = np.array(positions_single)
                s_single.fillna()
                s_single.update_direct()
                s_single_copy = copy.deepcopy(s_single)
                self._structures.append(s_single_copy)
                del s_single
                s_single = Structure()
                if(verbose):
                    print(f'\r{len(self._structures)}', end = '')
        if(verbose):
            print()
    
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
    
    def to_n2p2(self, file_name = "input.data", verbose = 1):
        with open(file_name, 'w') as f:
            if(verbose):
                print('Writing the n2p2 file...')
            
            num_structures = len(self._structures)
            for i in range(num_structures):
                f.write('begin\n')
                f.write('comment Created_by_datana\n')
                f.write(f'lattice {self._structures[i]._lattice[0][0]:.9f} {self._structures[i]._lattice[0][1]:.9f} {self._structures[i]._lattice[0][2]:.9f}\n')
                f.write(f'lattice {self._structures[i]._lattice[1][0]:.9f} {self._structures[i]._lattice[1][1]:.9f} {self._structures[i]._lattice[1][2]:.9f}\n')
                f.write(f'lattice {self._structures[i]._lattice[2][0]:.9f} {self._structures[i]._lattice[2][1]:.9f} {self._structures[i]._lattice[2][2]:.9f}\n')
                for j in range(len(self._structures[i]._species)):
                    f.write(f'atom {self._structures[i]._coords_cartesian[j][0]:.5f} {self._structures[i]._coords_cartesian[j][1]:.5f} {self._structures[i]._coords_cartesian[j][2]:.5f} ')
                    f.write(f'{self._structures[i]._species[j]} 0.0 0.0 ')
                    f.write(f'{self._structures[i]._forces[j][0]:.5f} {self._structures[i]._forces[j][1]:.5f} {self._structures[i]._forces[j][2]:.5f}\n')
                f.write(f'energy {self._structures[i]._energy:.8f}\n')
                f.write('charge 0.00\n')
                f.write('end\n')
                if(verbose):
                    print(f'\r{i+1}/{num_structures}', end = '')
        if(verbose):
            print()
    
    def to_POSCAR(self, verbose = 1, mode = 'Direct'):
        if(verbose):
            print('Writing the POSCAR file...')
        num_structures = len(self._structures)
        for i in range(num_structures):
            self._structures[i].to_POSCAR('POSCAR_'+str(i), mode = mode)
            if(verbose):
                print(f'\r{i+1}/{num_structures}', end = '')
        if(verbose):
            print()
    
    def to_xyz(self, file_name = 'data.xyz', verbose = 1):
        with open(file_name, 'w') as f:
            if(verbose):
                print('Writing the xyz file...')
            num_structures = len(self._structures)
            for i in range(num_structures):
                f.write(f'{len(self._structures[i]._species)}\n')
                f.write(f'Lattice=" {self._structures[i]._lattice[0][0]:.9f} {self._structures[i]._lattice[0][1]:.9f} {self._structures[i]._lattice[0][2]:.9f}')
                f.write(f' {self._structures[i]._lattice[1][0]:.9f} {self._structures[i]._lattice[1][1]:.9f} {self._structures[i]._lattice[1][2]:.9f}')
                f.write(f' {self._structures[i]._lattice[2][0]:.9f} {self._structures[i]._lattice[2][1]:.9f} {self._structures[i]._lattice[2][2]:.9f}"\n')
                for j in range(len(self._structures[i]._species)):
                    f.write(f'{self._structures[i]._species[j]}  {self._structures[i]._coords_cartesian[j][0]:.5f}  {self._structures[i]._coords_cartesian[j][1]:.5f}  {self._structures[i]._coords_cartesian[j][2]:.5f}\n')
                if(verbose):
                    print(f'\r{i+1}/{num_structures}', end = '')
        if(verbose):
            print()
