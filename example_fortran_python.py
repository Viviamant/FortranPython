from math import sqrt
from FortranPython import wrapfortran 
from numpy import random
import numpy as np 

myFortranLibrary = wrapfortran.FortranLibrary()

print("Test for Jacobi Diagonalization of a square symmetric Matrix")
A = [[1., sqrt(2.), 2.],
     [sqrt(2), 3., sqrt(2)],
     [2., sqrt(2.), 1.]]
myFortranLibrary.JacobiDiagonalization(A)
print()

print("Test on Least Square fitting of 1D x, y straight line")
rand = random.randn(6000)
X = [1.*i for i in range(6000)]
Y = [(2.*X[i] + rand[i]) for i in range(6000)]
myFortranLibrary.LeastSquareFit(X, Y)
print()

print("Test on Gauss Jordan Elimination of linear algebric equations solution")
A = [
     [ 2., 1.,-1., 2.],
     [ 4., 5.,-3., 6.],
     [-2., 5.,-2., 6.],
     [ 4.,11.,-4., 8.]]
C = [[ 5.,10.],
     [ 9.,18.],
     [ 4., 8.],
     [ 2., 4.]]
myFortranLibrary.GaussJordan(A, C)
print()

print("Test on Determinant calculation")
A = [ 
     [10., .0, -3.],
     [-2.,-4., 1.],
     [3., .0, 2.]
]
myFortranLibrary.GetDeterminant(A)
print()

print("Test on Fast Fourier Transform of complex array")
N = 1024 
T = 1. / 720.
x = np.linspace(.0, N*T, N)
f1 = 50.*2.*np.pi
f2 = 80.*2.*np.pi
signal = np.sin(f1*x) + .5*np.sin(f2*x) + 0j*x
myFortranLibrary.FastFourierTransform(x, signal)
print()

print("Test on Callback of a python function inside Fortran")
def getFunction(inputReal):
     x = inputReal
     return x*np.cos(10.*x**2)/(x**2 + 1.)
myFortranLibrary.Integrate(getFunction, .0, np.pi, 1000)
