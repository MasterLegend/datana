# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import cv2

class Figure:
    def __init__(self):
        self.flag = 1
    
    def __dfs_nb(self, data_binary, r, c):
        if(r == 0 or r == self.n_row - 1 or c == 0 or c == self.n_col - 1):
            self.flag = 0
            return
        if(data_binary[r][c] == False):
            return
        data_binary[r][c] = False
        self.__dfs_nb(data_binary, r - 1, c)
        self.__dfs_nb(data_binary, r + 1, c)
        self.__dfs_nb(data_binary, r, c - 1)
        self.__dfs_nb(data_binary, r, c + 1)
        return
    
    def __dfs_wb(self, data_binary, r, c):
        if(r < 0 or r > self.n_row - 1 or c < 0 or c > self.n_col - 1):
            return
        if(data_binary[r][c] == False):
            return
        data_binary[r][c] = False
        self.__dfs_wb(data_binary, r - 1, c)
        self.__dfs_wb(data_binary, r + 1, c)
        self.__dfs_wb(data_binary, r, c - 1)
        self.__dfs_wb(data_binary, r, c + 1)
        return
    
    def set_data(self, data):
        self.data = data
        self.n_row = len(self.data)
        self.n_col = len(self.data[0])
        return
    
    def get_data(self):
        return self.data
    
    def imread(self, path, mode = cv2.IMREAD_GRAYSCALE):
        self.data = cv2.imread(path, mode)
        self.n_row = len(self.data)
        self.n_col = len(self.data[0])
    
    def imshow(self):
        plt.imshow(self.data, 'gray')
        plt.xticks([])
        plt.yticks([])
        plt.show()
        
    def blur(self, ksize  = (1, 1)):
        self.data = cv2.blur(self.data, ksize)
        
    def count(self, threshold = 0.5, mode = 'avg'):
        data_binary = self.data > threshold * 255
        plt.imshow(data_binary, 'gray')
        res_nb = 0
        for r in range(self.n_row):
            for c in range(self.n_col):
                if(data_binary[r][c] == True):
                    self.flag = 1
                    self.__dfs_nb(data_binary, r, c)
                    res_nb += self.flag
        data_binary = self.data > threshold * 255
        res_wb = 0
        for r in range(self.n_row):
            for c in range(self.n_col):
                if(data_binary[r][c] == True):
                    self.flag = 1
                    self.__dfs_wb(data_binary, r, c)
                    res_wb += self.flag
        res = 0.5 * (res_nb + res_wb)
        print(f'The number of particles is {res}')
        return res
        
    def flatten(self, num_r = 3, num_c = 3):
        list_r = [int(i * self.n_row / num_r) for i in range(num_r + 1)]
        list_c = [int(i * self.n_col / num_c) for i in range(num_c + 1)]
        
        list_min = []
        for i in range(num_r * num_c):
            data_split = self.data[list_r[i//num_c]:list_r[i//num_c+1], list_c[i%num_c]:list_c[i%num_c+1]]
            index_min = np.argmin(data_split)
            r = int(index_min / len(data_split[0]))
            c = index_min % len(data_split[0])
            list_min.append([r + list_r[i//num_c], c + list_c[i%num_c], data_split[r][c]])
        list_min = np.array(list_min)
        
        def func_plane(coor, a, b, c):
            z = a*coor[0] + b*coor[1] + c
            return z
        coor = [list_min[:,0], list_min[:,1]]
        pfit, pcov = curve_fit(func_plane, coor, list_min[:,2])
        
        plane = []
        for r in range(self.n_row):
            for c in range(self.n_col):
                plane.append(pfit[0] * r + pfit[1] * c + pfit[2])
        plane = np.array(plane)
        plane = plane.reshape(self.n_row, self.n_col)
        
        data_flatten = self.data - plane
        data_min = np.min(data_flatten)
        data_max = np.max(data_flatten)
        minmax = (data_flatten - data_min) / (data_max - data_min)
        if(data_min < 0): data_min = 0
        if(data_max > 255): data_max = 255
        data_flatten = minmax * (data_max - data_min) + data_min
        
        fig_flatten = Figure()
        fig_flatten.set_data(data_flatten)
        
        return fig_flatten
    
    def flatten3p(self, list_coor):  
        x1 = list_coor[0][0]
        y1 = list_coor[0][1]
        z1 = int(self.data[x1][y1])
        x2 = list_coor[1][0]
        y2 = list_coor[1][1]
        z2 = int(self.data[x2][y2])
        x3 = list_coor[2][0]
        y3 = list_coor[2][1]
        z3 = int(self.data[x3][y3])
        
        para_a = (y2 - y1)*(z3 - z1) - (z2 -z1)*(y3 - y1)
        para_b = (x3 - x1)*(z2 - z1) - (x2 - x1)*(z3 - z1)
        para_c = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)
        para_d = -(para_a * x1 + para_b * y1 + para_c * z1)
        
        plane = []
        for r in range(self.n_row):
            for c in range(self.n_col):
                plane.append(-(para_a * r + para_b * c + para_d) / para_c)
        plane = np.array(plane)
        plane = plane.reshape(self.n_row, self.n_col)
        
        data_flatten = self.data - plane
        data_min = np.min(data_flatten)
        data_max = np.max(data_flatten)
        minmax = (data_flatten - data_min) / (data_max - data_min)
        if(data_min < 0): data_min = 0
        if(data_max > 255): data_max = 255
        data_flatten = minmax * (data_max - data_min) + data_min
        
        fig_flatten = Figure()
        fig_flatten.set_data(data_flatten)
        
        return fig_flatten