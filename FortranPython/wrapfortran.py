import ctypes as ct
import numpy as np
import sys
from time import perf_counter 
from .utils import pointer_cmplx_array, pointer_array

class FortranLibrary:
	def __init__(self, fort_so_path='/usr/local/lib/mylib.so'):
		self.fort_so_path = fort_so_path
		try:
			self.fort_lib = ct.CDLL(self.fort_so_path)
			print(f"Fortran library object '{self.fort_so_path}' is loaded in memory")
			print()
		except:
			print(f'Exception in Fortran shared object {self.fort_so_path}; ', end='')
			print('Maybe a compilation error has occurred or the file could not be found')
		sys.stdout.flush()

	def Integrate(self, f, r0, r1, i):
		Integrate_ = self.fort_lib.Integrate
		getFunVal_proc = ct.CFUNCTYPE(ct.c_double, ct.c_double)
		getFoo_pntr = getFunVal_proc(f)
		x0 = ct.c_double(r0)
		x1 = ct.c_double(r1)
		outputReal = ct.c_double(.0)
		ii = ct.c_int(i)
		start = perf_counter()
		Integrate_(getFoo_pntr, ct.byref(x0), ct.byref(x1), ct.byref(ii), ct.byref(outputReal))
		end = perf_counter()
		print(f"value of outputReal received in Python for {ii.value} intervals: {np.double(outputReal):.7f}")
		print(f'Elapsed time in Fortran {(end - start)*1e3} ms')
		sys.stdout.flush()

	def GetDeterminant(self, A):
		GetDeterminant_ = self.fort_lib.GetDeterminant 
		GetDeterminant_.restype = None 
		a = np.array(A, order='F')
		shape = a.shape
		assert shape[0] == shape[1], 'A must be a Square Matrix'
		n = ct.c_int(shape[0])
		det = ct.c_double()
		start = perf_counter()
		GetDeterminant_(ct.byref(n), pointer_array(a), ct.byref(det))
		end = perf_counter()
		print(f'Determinant of matrix \n{np.array(A)} \nis {det.value}')
		print(f'Elapsed time in Fortran {(end - start)*1e3} ms')
		sys.stdout.flush()

	def FastFourierTransform(self, x, S, plotting=False):
		FastFourierTransform_ = self.fort_lib.FastFourierTransform
		FastFourierTransform_.restype = None 
		f = np.array(S, order='F')
		n = ct.c_int(len(S))
		print('Original Signal:')
		print(np.array(S))
		start = perf_counter()
		FastFourierTransform_(	ct.byref(n), 
								pointer_cmplx_array(f))
		end = perf_counter()
		print('Fourier Transform:')
		print(f)
		if plotting:
			import matplotlib.pyplot as plt 
			print('Plotting...')
			T = x[1] - x[0]    
			N = len(S)
			freq = np.linspace(.0, 1.//T, N)
			_, axes = plt.subplots(2)
			axes[0].plot(x, S.real)
			axes[0].plot(x, S.imag)
			axes[1].plot(freq[:N//2], 2./N * abs(f[:N//2]))
			plt.show()
		print(f'Elapsed time in Fortran {(end - start)*1e3} ms')
		sys.stdout.flush()

	def JacobiDiagonalization(self, A):
		JacobiDiagonalization_ = self.fort_lib.JacobiDiagonalization 
		JacobiDiagonalization_.restype = None 
		shape = np.array(A).shape
		assert shape[0] == shape[1], 'Matrix A must be a square 2d array'
		a = np.array(A, order='F')
		q = np.zeros(shape, order='F')
		w = np.zeros(shape[0], order='F')
		start = perf_counter()
		JacobiDiagonalization_( pointer_array(a), 
								pointer_array(q), 
								pointer_array(w))
		end = perf_counter()
		print('Original Matrix:')
		print(np.array(A))
		print('Eigensystem:')
		eigenvectors = q.T
		eigenvalues = w
		for eigenvalue, eigenvector in zip(eigenvalues, eigenvectors):
			print(eigenvector, eigenvalue)
		print(f'Elapsed time in Fortran {(end - start)*1e3} ms')
		sys.stdout.flush()

	def LeastSquareFit(self, X, Y):
		LeastSquareFit_ = self.fort_lib.LeastSquareFit
		LeastSquareFit_.restype = None 
		x = np.array(X, order='F')
		y = np.array(Y, order='F')
		n1 = x.shape[0]
		n2 = y.shape[0]
		assert n1 == n2, 'X and Y must have the same length'
		n = ct.c_int(n1)
		m = ct.c_double()
		q = ct.c_double()
		erm = ct.c_double()
		erq = ct.c_double()
		r2 = ct.c_double()
		start = perf_counter()
		LeastSquareFit_(ct.byref(n), 
						pointer_array(x), 
						pointer_array(y), 
						ct.byref(m), 
						ct.byref(q), 
						ct.byref(r2), 
						ct.byref(erm), 
						ct.byref(erq))
		end = perf_counter()
		M = m.value 
		ErM = erm.value 
		Q = q.value 
		ErQ = erq.value 
		R2 = r2.value 
		print(f'Angular coefficient= {M:e}+/-{ErM:e}')
		print(f'Intercept= {Q:e}+/-{ErQ:e}')
		print(f'r^2= {R2:e}')
		print(f'Elapsed time in Fortran {(end - start)*1e3} ms')
		sys.stdout.flush()

	def GaussJordan(self, A, B):
		GaussJordan_ = self.fort_lib.GaussJordan 
		GaussJordan_.restype = None
		a = np.array(A, order='F')
		b = np.array(B, order='F')
		na = a.shape
		assert na[0] == na[1], 'A must be a square matrix'
		nb = b.shape
		assert na[1] == nb[0], 'A 2nd dimension and B 1st dimension must be equal'
		n = ct.c_int(nb[0]) 
		m = ct.c_int(nb[1])
		start = perf_counter()
		GaussJordan_(ct.byref(n), 
					 ct.byref(m), 
					 pointer_array(a), 
					 pointer_array(b))
		end = perf_counter()
		invA = a 
		sol = b
		print('Coefficients matrix:')
		print(np.array(A))
		print('RHS column vectors:')
		print(np.array(B))
		print('Solution column vectors:')
		print(sol)
		print('Inverse coefficients Matrix:')
		print(invA)
		print(f'Elapsed time in Fortran {(end - start)*1e3} ms')
		sys.stdout.flush()
