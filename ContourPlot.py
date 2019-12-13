#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 13:26:25 2019

@author: molly199
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

import pandas as pd

#Import the .plt file
test = pd.read_csv('UVP.plt',skiprows=2,delimiter='\t')
#Convert to numpy array
test=test.values

#Set the variables as arrays
x = test[:,0]
y = test[:,1]
u = test[:,2]
v = test[:,3]
p = test[:,4]

#Contour Plot (taken from Stack Overflow and edited to fit output)
def plot_contour(x,y,z,resolution = 100, contour_method = 'linear'):
    resolution = str(resolution)+'j'
    X,Y = np.mgrid[min(x):max(x):complex(resolution), min(y):max(y):complex(resolution)]
    points = [[a,b] for a,b in zip(x,y)]
    Z = griddata(points,z,(X,Y), method=contour_method)
    return X,Y,Z

#Use contour plot function to convert data to 2d array
Xp,Yp,P = plot_contour(x,y,p)

#Plot Pressure
plt.contourf(Xp,Yp,P)
plt.title('Contour Plot of Pressure')
plt.savefig('PressureContour.png')
plt.show()

#Use contour plot function to convert data 2d array
Xu,Yu,U = plot_contour(x,y,u)

#Plot U
plt.contourf(Xu,Yu,U)
plt.title('Contour Plot of U Velocity')
plt.savefig('VelocityContour.png')
plt.show()


