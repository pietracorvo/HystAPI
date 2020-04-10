#!/usr/bin/env python

import HystAPI
import matplotlib.pyplot as plt

# Load five data sets from the folder /testdata 
sample1 = HystAPI.Curve('./testData/FeSi3_annealed', ' Fesi3 (annealed)')
sample2 = HystAPI.Curve('./testData/FeSi3_RT', 'FeSi3 (T = 77K)')
sample3 = HystAPI.Curve('./testData/FeSi3_RT_ht', 'FeSi3 (T = 77K, heat treated)')
sample4 = HystAPI.Curve('./testData/FeSi3_LN', ' FeSi3 (T = 293K)')
sample5 = HystAPI.Curve('./testData/FeSi3_LN_ht', 'FeSi3 (T = 293K, heat treated)')

# Initialize a measurement ...
measurement = HystAPI.Measurement('FeSi3 (Electrical Steel)')

# ... and add datasets to the measurement 
measurement.addCurve(sample1)
measurement.addCurve(sample2)
measurement.addCurve(sample3)
measurement.addCurve(sample4)
measurement.addCurve(sample5)

# Plot the hystereses of each sample
sample1.plotHyst()
sample2.plotHyst()
sample3.plotHyst()
sample4.plotHyst()
sample5.plotHyst()

# Plot the main results from all samples together
measurement.plotLosses()
measurement.plotHc()
measurement.plotMuMax()

# Save calculated fit parameters in a text file 
measurement.printInfo('testFitParameters.txt')

plt.show()

