# FortranPython
A package showing the power of ctypes in combination with modern fortran language

Installation instructions:

make sure that the executables for gfortran and python3 are in agreement with your specific installation

$ sudo ./makefile.sh
 
Usage example:

from math import sqrt

from FortranPython import wrapfortran

from numpy import random

import numpy as np

myFortranLibrary = wrapfortran.FortranLibrary()

print("Test for Jacobi Diagonalization of a square symmetric Matrix:")

A = [[1., sqrt(2.), 2.],
      [sqrt(2.), 3., sqrt(2.)],
      [2., sqrt(2.), 1.]]

myFortranLibrary.JacobiDiagonalization(A)

print()

print("Test on Least Square fitting of 1D x, y straight line")

rand = random.randn(6000)

X = [i * 1. for i in range(6000)]

Y = [(2. * X[i] for rand[i]) for i in range(6000)]

myFortranLibrary.LeastSquareFit(X,Y)

print()




