#!/usr/bin/env python
# coding: utf8

###############################################################################

# DESCRIPTION: Module to analyse and plot Hysteresis data
#              The ojects Hyst, Curve and Measurement represent the data of
#              a single hysteresis object (Hyst), a collection of multiple
#              hystereses measured at different frequencies (Curve) and a
#              collection of multiple Curve objects for different samples
#              (Measurement)
# AUTHOR       Rabensteiner Alexander
# DATE:        04.04.2019

###############################################################################

import os, sys

import matplotlib
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt

import numpy as np
from scipy.signal import savgol_filter as smoothing
from scipy.optimize import curve_fit

from prettytable import PrettyTable

class Hyst():
    '''Contains data of one hysteresis measurement'''

    # Constructor
    def __init__(self, pathToData='', freq=-1.0):
        # Initialize all members
        self.pathToData = pathToData
        self.freq = freq
        self.extrema = []   # array indices of minima and maxima in H (has 4 entries in this order: [max1, min1, max2, min2])
        self.data = {'t': [], 'H': [], 'J': []}  # container for the raw data
        self.muDiff, self.muAbs = ([] for i in range(2))   # arrays of calculated data
        self.losses, self.Hc, self.Jr, self.muMax = (0.0 for i in range(4)) # scalars of calculated data
        # Load data file from pathToData in the class
        self._setData()
        # Calculate secondary data (losses, Hc, etc.) from raw data
        self._setCalc()

    # Set data
    def _setData(self):
        # Open filestream with data
        try:
            fileStream = open(self.pathToData, 'r')
        except FileNotFoundError:
            print('ERROR: The file '+self.pathToData+' does not exist!')
            sys.exit(1)
        # Parse header (containing information like frequency, etc.) and hysteresis data seperately
        headerLineNum = 27  # that is a hard coded number, the number of header lines
        header = self._readHeader(fileStream, headerLineNum) # the header contains information like the frequency
        # If frequency is not given to constructor parse it from header of data file
        if self.freq==-1.0:
            freqDummy = header[15][4:-1] # for this specific files the frequency is in line 15
            try:
                self.freq = float(freqDummy)
            except ValueError:
                print('ERROR: The parsed frequency from file '+self.pathToData+' is not a float!')
                sys.exit(1)
        # Read time, applied field and magnetic polarisation from data file
        t, H, J = self._readData(fileStream, headerLineNum)
        self.data = { 't': t, 'H': H, 'J': J }
        # Only the hysteresis between first and second maximum is interesting, that means we need the extrema of H
        extrema = self._findExtrema(self.data['H'])
        self.extrema = [extrema[0][0], extrema[1][0], extrema[0][1], extrema[1][1]] # [max1, min1, max2, min2]
    def _readHeader(self, fileStream, headerLineNum):
        header = []
        for i in range(headerLineNum):
            header.append(fileStream.readline())
        return header
    def _readData(self, fileStream, headerLineNum):
        raw_data = fileStream.readlines()
        new_data = [ el.split('\t') for el in raw_data[28:] ]
        t = [ float(line[0]) for line in new_data ]
        H = [ float(line[3]) for line in new_data ]
        J = [ float(line[9]) for line in new_data ]
        return np.asarray(t), np.asarray(H), np.asarray(J)
    def _setCalc(self):
        self.losses = self._calcLosses()
        self.muDiff = self._calcMuDiff()
        self.muAbs = self._calcMuAbs()
        self.muMax = np.max(self.muAbs)
        self.Hc = self._calcHc()
        self.Jr = self._calcJr()

    # Utilities
    def _calcLosses(self, steps=20):
        # Get losses by surface under hysteresis (first smooth data, then integrate by trapezoidrule)
        maxima, minima = self._findExtrema(self.data['H'])
        data_x = smoothing(np.copy(self.data['H'])-np.min(self.data['H']), 501, 4)
        data_y = smoothing(np.copy(self.data['J'])-np.min(self.data['J']), 501, 4)
        neg_val = np.abs(np.trapz(data_y[maxima[0]:minima[0]], data_x[maxima[0]:minima[0]], steps))
        pos_val = np.abs(np.trapz(data_y[minima[0]:maxima[1]], data_x[minima[0]:maxima[1]], steps))
        return np.abs(pos_val-neg_val)
    def _calcMuDiff(self):
        # Differentail permeability (is very noise sensitive)
        maximum = self.extrema[0]
        dy = [ np.abs(np.min(self.data['J'][:int(0.5*maximum)])) if np.min(self.data['J'][:int(0.5*maximum)])<0 else 0 ]
        dx = [ np.abs(np.min(self.data['H'][:int(0.5*maximum)])) if np.min(self.data['H'][:int(0.5*maximum)])<0 else 0 ]
        Jsm = smoothing(dy + self.data['J'][:maximum], 501, 3)
        Hsm = smoothing(dx + self.data['H'][:maximum], 501, 3)
        mu_diff = np.empty(0)
        mu_H = np.empty(0)
        for i in range(Hsm.size-1):
            if Hsm[i]<50:
                mu = 0
            else:
                mu = (Jsm[i+1] - Jsm[i]) / ((Hsm[i+1]-Hsm[i])*1.2566E-6)
            mu_diff = np.append(mu_diff, mu)
        mu_diff[mu_diff<0] = 0
        return mu_diff
    def _calcMuAbs(self,):
        # Absolute permeability: Slope of straight line connecting origin and value pair (H,J)
        maximum = self.extrema[0]
        mu_abs = np.empty(0)
        dy = [ np.abs(np.min(self.data['J'][:int(0.5*maximum)])) if np.min(self.data['J'][:int(0.5*maximum)])<0 else 0 ]
        dx = [ np.abs(np.min(self.data['H'][:int(0.5*maximum)])) if np.min(self.data['H'][:int(0.5*maximum)])<0 else 0 ]
        Jsm = smoothing(dy + self.data['J'][:maximum], 501, 3)
        Hsm = smoothing(dx + self.data['H'][:maximum], 501, 3)
        for i in range(Hsm.size):
            if Hsm[i]<50:
                mu = 0
            else:
                mu = Jsm[i]/(Hsm[i]*1.2566E-6)
            mu_abs = np.append(mu_abs, mu)
        mu_abs[mu_abs<0] = 0
        return mu_abs
    def _calcHc(self):
        # Coercitivity field [A/m]
        maxima, minima = self._findExtrema(self.data['H'])
        Hc1 = self.data['H'][maxima[0]+np.argmin(np.abs(self.data['J'][maxima[0]:minima[0]]))]
        Hc2 = self.data['H'][minima[0]+np.argmin(np.abs(self.data['J'][minima[0]:maxima[1]]))]
        Hc = 0.5*(np.abs(Hc1)+np.abs(Hc1))
        return Hc
    def _calcJr(self):
        # Remanencepolarisation [T]
        maxima, minima = self._findExtrema(self.data['H'])
        Jr1 = self.data['J'][maxima[0]+np.argmin(np.abs(self.data['H'][maxima[0]:minima[0]]))]
        Jr2 = self.data['J'][minima[0]+np.argmin(np.abs(self.data['H'][minima[0]:maxima[1]]))]
        Jr = 0.5*(np.abs(Jr1)+np.abs(Jr2))
        return Jr
    def _findExtrema(self, data):
        offset =  int(0.01*len(data)) # the default value is problematic for small number of points
        ind_max = np.empty(0, dtype='intp')
        ind_min = np.empty(0, dtype='intp')
        maximum = -1000000000
        minimum = 1000000000
        flag = 1
        i = offset
        while i<data.size:
            if flag==1:
                tmp_max = np.max(data[i-offset:i+offset])
                if tmp_max>=maximum:
                    maximum = tmp_max
                    dummy_ind_max = np.argmax(data[i-offset:i+offset])+ i - offset
                elif tmp_max<maximum:
                    maximum = -1000000000
                    ind_max = np.append(ind_max, int(dummy_ind_max))
                    flag = -1
            elif flag==-1:
                tmp_min = np.min(data[i-offset:i+offset])
                if tmp_min<=minimum:
                    minimum = tmp_min
                    dummy_ind_min = np.argmin(data[i-offset:i+offset])+ i - offset
                elif tmp_min>minimum:
                    minimum = 1000000000
                    ind_min = np.append(ind_min, int(dummy_ind_min))
                    flag = +1
            i += offset
            ind_max = ind_max[0:2]
            ind_min = ind_min[0:2]
        return ind_max, ind_min

    # Output
    def plotHyst(self, title=''):
        # Build template axes
        bigLine = 1.5
        fig = plt.figure(figsize=(10,7))
        ax = fig.add_subplot(111)
        ax.set_title(title)
        ax.set_ylabel('J [T]')
        ax.set_ylabel('H [A/m]')
        ax.set_ylim(-2, 2)
        ax.set_xlim(-10000, 10000)
        ax.spines['top'].set_linewidth(bigLine)
        ax.spines['bottom'].set_linewidth(bigLine)
        ax.spines['left'].set_linewidth(bigLine)
        ax.spines['right'].set_linewidth(bigLine)
        ax.tick_params(which='major', length=8, direction='in', width=bigLine, pad=10, top='on', right='on')
        ax.set_yticks(np.linspace(-2,2,11))
        ax.set_xticks(np.linspace(-10000,10000,11))
        ax.tick_params(which='minor', length=4, direction='in', width=bigLine, pad=10, top='on', right='on')
        ax.set_yticks(np.linspace(-2,2,21), minor=True)
        ax.set_xticks(np.linspace(-10000,10000,51), minor=True)
        ax.grid(which='major')
        ax.grid(which='minor', color='xkcd:light grey', linestyle=':')
        ax.axhline(linewidth=bigLine, color='black')
        ax.axvline(linewidth=bigLine, color='black')
        ax = fig.axes[0]
        # And add data to the plot
        ax.plot(self.data['H'], self.data['J'], label=str(int(self.freq))+' Hz')
        ax.legend(loc=4)
        return fig
    def showInfo(self):
        print('Hysteresis information:')
        print('----------------------')
        print('Data loaded from: '+self.pathToData)
        print('Frequency [Hz]: %.5f' % self.freq)
        print('Losses [W/m^3]: %.5f' % self.losses)
        print('Hc [A/m]: %.5f' % self.Hc)
        print('Jr [T]: %.5f' % self.Jr)
        print('MuMax []: %.5f' % self.muMax)
        print('\n')

class Curve():
    '''A Curve object collects the data of the same sample measured with different frequencies'''

    # Constructor
    def __init__(self, dataLocation='', curveName=''):
        self.dataLocation = dataLocation # this is a folder containing data files of measurement at different frequencies for the same sample
        self.curveName = curveName # name used default as title of plots
        self.data = {}
        # Load the data of any hysteresis inside the folder dataLocation
        self.data = self._setData(dataLocation)
        # Calculate fits from given data
        self._updateData()

    # Inject Hyst objects into Curve (also an empty Curve() can be build an Hyst() objects injected by hand)
    def addHyst(self, hyst):
        self.data[hyst.freq] = hyst
        self._updateData()

    # Set data
    def _setData(self, basePath):
        # Get the names of all files in basePath
        fileList = self._getFilenames(basePath)
        data = {}
        # Write Hyst instances into a list
        for fileName in fileList:
            dummyInstance = Hyst(fileName)
            data[dummyInstance.freq] = dummyInstance
        return data
    def _getFilenames(self, basePath):
        filePaths = []
        # Check if path exists
        if os.path.isdir(basePath)==false:
            print('ERROR: '+basePath+' is not a path!')
            sys.exit(1)
        # Pathwalk with the given directory to get filenames in it
        for dirpath, dirnames, filenames in os.walk(basePath):
            for filename in filenames:
                filePaths.append(os.path.join(dirpath, filename))
        return filePaths
    def _updateData(self):
        # Calculate the fits for data
        self.lossesFit = self._calcFit( self._freqs(), self._losses() )
        self.HcFit = self._calcFit( self._freqs(), self._Hc() )
        self.JrFit = self._calcFit( self._freqs(), self._Jr() )
        self.muMaxFit = self._calcFitExp( self._freqs(), self._muMax() )
    def _calcFit(self, freq, y):
        # Calculate generic fit (for losses, Hc, etc.)
        def fitfun(x, a, b, c):
            return a+b*x+c*np.sqrt(x)
        y = [x for _,x in sorted(zip(freq, y))]
        freq = sorted(freq)
        param, dummy = curve_fit(fitfun, freq, y)
        evaluation = fitfun(np.linspace(np.min(freq), np.max(freq), 1000), *param)
        xfit = np.linspace(np.min(freq), np.max(freq), 1000)
        return { 'val':y, 'fit':evaluation, 'x':xfit, 'par':param }
    def _calcFitExp(self, freq, y):
        # calculate exponential fit (e.g. for permeability muAbs)
        def fitfun(x, y0, A1, A2, t1, t2):
            return y0 + A1 * np.exp(-x/t1) + A2 * np.exp(-x/t2)
        y = [x for _,x in sorted(zip(freq, y))]
        freq = sorted(freq)
        param, dummy = curve_fit(fitfun, freq, y)
        evaluation = fitfun(np.linspace(np.min(freq), np.max(freq), 1000), *param)
        xfit = np.linspace(np.min(freq), np.max(freq), 1000)
        return { 'val':y, 'fit':evaluation, 'x':xfit, 'par':param }

    # Getter
    def _freqs(self):
        return sorted(self.data.keys())
    def _losses(self):
        return [ self.data[key].losses for key in self._freqs() ]
    def _Hc(self):
        return [ self.data[key].Hc for key in self._freqs() ]
    def _Jr(self):
        return [ self.data[key].Jr for key in self._freqs() ]
    def _muMax(self):
        return [ self.data[key].muMax for key in self._freqs() ]

    # Output
    def plotHyst(self):
        # Build template axes
        # Build template axes
        bigLine = 1.5
        fig = plt.figure(figsize=(10,7))
        ax = fig.add_subplot(111)
        ax.set_title(self.curveName)
        ax.set_ylabel('J [T]')
        ax.set_ylabel('H [A/m]')
        ax.set_ylim(-2, 2)
        ax.set_xlim(-10000, 10000)
        ax.spines['top'].set_linewidth(bigLine)
        ax.spines['bottom'].set_linewidth(bigLine)
        ax.spines['left'].set_linewidth(bigLine)
        ax.spines['right'].set_linewidth(bigLine)
        ax.tick_params(which='major', length=8, direction='in', width=bigLine, pad=10, top='on', right='on')
        ax.set_yticks(np.linspace(-2,2,11))
        ax.set_xticks(np.linspace(-10000,10000,11))
        ax.tick_params(which='minor', length=4, direction='in', width=bigLine, pad=10, top='on', right='on')
        ax.set_yticks(np.linspace(-2,2,21), minor=True)
        ax.set_xticks(np.linspace(-10000,10000,51), minor=True)
        ax.grid(which='major')
        ax.grid(which='minor', color='xkcd:light grey', linestyle=':')
        ax.axhline(linewidth=bigLine, color='black')
        ax.axvline(linewidth=bigLine, color='black')
        ax = fig.axes[0]
        # And add data to the plot
        for freq in self._freqs():
            ax.plot(self.data[freq].data['H'], self.data[freq].data['J'], label=str(int(freq))+' Hz')
        ax.legend(loc=4)
        return fig
    def _plotTemplate(self, xlabel='', ylabel='', title='', data=[], fitx=[], fity=[]):
        # wrapper for the production of figures
        fig = plt.figure(figsize=(10,7))
        ax = fig.add_subplot(111)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(color='0.75', linestyle='--', linewidth=0.5)
        ax = fig.axes[0]
        ax.plot(self._freqs(), data, marker='o', linestyle='', color='blue')
        ax.plot(fitx, fity, linestyle='--', color='blue')
        return fig
    def plotLosses(self):
        fig = self._plotTemplate('f [1/s]', 'P/V [W/m$^3$]', 'Loss curve '+self.curveName,
                            self._losses(), self.lossesFit['x'], self.lossesFit['fit'])
        return fig
    def plotHc(self):
        fig = self._plotTemplate('f [1/s]', 'H [A/m]', 'Coercive Field H$_{c}$ '+self.curveName,
                            self._Hc(), self.HcFit['x'], self.HcFit['fit'])
        return fig
    def plotMuMax(self):
        fig = self._plotTemplate('f [1/s]', '$\mu_{max}$ [H/m]', 'Maximal permeability $\mu_{max}$ '+self.curveName,
                            self._muMax(), self.muMaxFit['x'], self.muMaxFit['fit'])
        return fig
    def showInfo(self):
        print('Fit parameters:')
        print('---------------')
        print('Losses: W(f) = %.3f + %.3f *  + %.3f * ^0.5' % (self.lossesFit['par'][0], self.lossesFit['par'][1], self.lossesFit['par'][2]))
        print('Coercive field Hc(f) = %.3f + %.3f * f + %.5f * f^0.5' % (self.HcFit['par'][0], self.HcFit['par'][1], self.HcFit['par'][2]))
        print('Maximal permeability: muMax(f) = %.3f + %.3f * e^(-f/%.3f) + %.3f * e^(-f/%.3f)' % (self.muMaxFit['par'][0], self.muMaxFit['par'][1], self.muMaxFit['par'][2], self.muMaxFit['par'][3], self.muMaxFit['par'][4]))

class Measurement():
    ''' The Measurement object collects multiple Curve objects to plot them all together'''

    # Constructor
    def __init__(self, measurementName=''): # the constructor always builds an empty instance which must be filled with the class method .addCurve()
        self.measurementName = measurementName
        self.collection = {}    # that is an empty container for multiple Curve objects

    # Inject Curve objects in Measurement
    def addCurve(self, curve):
        self.collection[curve.curveName] = curve

    # Output
    def _plotTemplate(self, xlabel='', ylabel='', title=''):
        # formatting
        fig = plt.figure(figsize=(10,7))
        ax = fig.add_subplot(111)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(color='0.75', linestyle='--', linewidth=0.5)
        return fig
    def plotLosses(self):
        fig = self._plotTemplate('f [1/s]', 'P/V [W/m$^3$]', 'Loss curve '+self.measurementName)
        ax = fig.axes[0]
        # And data to plot
        colormap = plt.cm.nipy_spectral(np.linspace(0.25, 0.95, len(self.collection.keys())))
        for curveName, col in zip(self.collection.keys(), colormap):
            ax.plot(self.collection[curveName]._freqs(), self.collection[curveName]._losses(), marker='o', linestyle='', label=curveName, color=col)
            ax.plot(self.collection[curveName].lossesFit['x'], self.collection[curveName].lossesFit['fit'], linestyle='--', color=col)
        ax.legend(loc=4)
        return fig
    def plotHc(self):
        fig = self._plotTemplate('f [1/s]', 'H [A/m]', 'Coercive Field H$_{c}$ '+self.measurementName)
        ax = fig.axes[0]
        # Add data to plot
        colormap = plt.cm.nipy_spectral(np.linspace(0.25, 0.95, len(self.collection.keys())))
        for curveName, col in zip(self.collection.keys(), colormap):
            ax.plot(self.collection[curveName]._freqs(), self.collection[curveName]._Hc(), marker='o', linestyle='', label=curveName, color=col)
            ax.plot(self.collection[curveName].HcFit['x'], self.collection[curveName].HcFit['fit'], linestyle='--', color=col)
        ax.legend(loc=4)
        return fig
    def plotMuMax(self):
        fig = self._plotTemplate('f [1/s]', '$\mu_{max}$ [H/m]', 'Maximal permeability $\mu_{max}$ '+self.measurementName)
        ax = fig.axes[0]
        # Add data to plot
        colormap = plt.cm.nipy_spectral(np.linspace(0.25, 0.95, len(self.collection.keys())))
        for curveName, col in zip(self.collection.keys(), colormap):
            ax.plot(self.collection[curveName]._freqs(), self.collection[curveName]._muMax(), marker='o', linestyle='', label=curveName, color=col)
            ax.plot(self.collection[curveName].muMaxFit['x'], self.collection[curveName].muMaxFit['fit'], linestyle='--', color=col)
        return fig
    def _buildTable(self):
        # Make the data table
        table = PrettyTable()
        table.field_names = ['Sample',
                             'W_stat [J/m^3]', 'W_cl [s*J/m^3]', 'W_ex [s^(1/2)*J/m^3]',
                             'Hc_stat [A/m]', 'Hc_cl [s*A/m]', 'Hc_ex [s^(1/2)*A/m]']
        # And plug values
        for sample in self.collection.keys():
            table.add_row([self.collection[sample].curveName,
                          self.collection[sample].lossesFit['par'][0],
                          self.collection[sample].lossesFit['par'][1],
                          self.collection[sample].lossesFit['par'][2],
                          self.collection[sample].HcFit['par'][0],
                          self.collection[sample].HcFit['par'][1],
                          self.collection[sample].HcFit['par'][2]])
        return table.get_string()
        # Write to terminal od to file if fileName is specified
    def printInfo(self, fileName=''):
        tableString = self._buildTable()
        # Select name of file for output
        outputFileName = self.measurementName+'fitParameters.txt'
        if fileName!='':
            outputFileName = fileName
        # Write table to file
        fileStream = open(outputFileName, 'w')
        fileStream.write(tableString)
        fileStream.close()
