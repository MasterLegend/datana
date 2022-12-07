# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

class Spectrum1D:
    def __init__(self, x_values, y_values):
        '''
        Constructor
        '''
        self._x_values = x_values
        self._y_values = y_values
    
    def __normal_distribution(self, x, mean, sigma):
        '''
        Normal distribution with mean as mean and sigma as sigma
        '''
        return np.exp(-1 * ((x - mean) ** 2) / (2 * (sigma ** 2))) / (np.sqrt(2 * np.pi) * sigma)
    
    def get_values(self):
        '''
        Get x and y values
        '''
        return self._x_values, self._y_values
    
    def set_values(self, x_values, y_values):
        '''
        Set the x and y values manually
        '''
        self._x_values = x_values
        self._y_values = y_values
    
    def smooth(self, method = 'sg', window_length = 9, polyorder = 2):
        y_sg = savgol_filter(self._y_values, window_length, polyorder)
        return Spectrum1D(self._x_values, y_sg)
    
    def find_peak(self, window_length = 5, verbose = 1):
        list_peak = []
        x_min = np.min(self._x_values)
        x_max = np.max(self._x_values)
        x_dx = (x_max - x_min) / len(self._x_values)
        window_step = int(window_length / x_dx)
        for i in range(len(self._y_values)):
            if(self._x_values[i] < x_min + window_length or self._x_values[i] > x_max - window_length):
                continue
            if(self._y_values[i] == np.max(self._y_values[i - window_step : i + window_step])):
                dic_peak = {'x':self._x_values[i], 'height':self._y_values[i]}
                list_peak.append(dic_peak)
        df_peak = pd.DataFrame(list_peak)
        df_peak = df_peak.sort_values('height', ascending = False, ignore_index = True)
        if(verbose):
            if(df_peak.empty):
                print('No peak!')
            else:
                print(f'{len(df_peak)} peak(s) is found')
        return df_peak
    
    def sub_baseline(self, method = 'linear', x_min = 0, x_max = 0):
        if(x_min == 0 and x_max == 0):
            x_min = np.min(self._x_values)
            x_max = np.max(self._x_values)
        dist_min = np.abs(self._x_values - x_min)
        dist_max = np.abs(self._x_values - x_max)
        index_min = np.argsort(dist_min)[0]
        index_max = np.argsort(dist_max)[0]
        slope = (self._y_values[index_max] - self._y_values[index_min]) / (self._x_values[index_max] - self._x_values[index_min])
        intercept = self._y_values[index_max] - slope * self._x_values[index_max]
        bl = slope * self._x_values + intercept
        y_sub = self._y_values - bl
        y_sub -= np.min(y_sub)
        return Spectrum1D(self._x_values, y_sub)
    
    def cal_integral(self, x_min = 0, x_max = 0, sub_bl = True):
        if(x_min == 0 and x_max == 0):
            x_min = self._x_values[0]
            x_max = self._x_values[len(self._x_values) - 1]
        dist_min = np.abs(self._x_values - x_min)
        dist_max = np.abs(self._x_values - x_max)
        index_min = np.argsort(dist_min)[0]
        index_max = np.argsort(dist_max)[0]
        area = 0.0
        for i in range(index_min, index_max):
            area += (self._y_values[i] + self._y_values[i + 1]) * (self._x_values[i + 1] - self._x_values[i]) * 0.5
        if(sub_bl):
            area -= (self._y_values[index_min] + self._y_values[index_max]) * (self._x_values[index_max] - self._x_values[index_min]) * 0.5
        print(f'area value is {area:.4f}')
        return area
        
    def peak_separation(self, num_peak = 1):
        para = np.ones(3 * num_peak)
        spec_peak = self.find_peak(verbose = 0)[:num_peak]
        for item in range(num_peak):
            para[3 * item] = spec_peak.height[item]
            para[3 * item + 1] = spec_peak.x[item]
        def func(x, *arg):
            res = 0
            for i in range(num_peak):
                res += arg[3*i] * self.__normal_distribution(x, arg[3*i+1], arg[3*i+2])
            return res
        pfit, pcov = curve_fit(func, self._x_values, self._y_values, p0 = para, bounds=(0, np.inf))

        res_fit = func(self._x_values, *pfit)
        plt.plot(self._x_values, self._y_values, 'k')
        plt.plot(self._x_values, res_fit, 'r--')
        for i in range(num_peak):
            plt.fill_between(self._x_values, 0, pfit[3*i] * self.__normal_distribution(self._x_values, pfit[3*i+1], pfit[3*i+2]), alpha = 0.5)
        plt.show()
        pfit_peak = [np.round(pfit[i], 4) for i in range(1, len(pfit), 3)]
        print(f'peak: {pfit_peak}')
        return
    
    def plot(self):
        plt.plot(self._x_values, self._y_values, 'k')
        plt.show()
        return