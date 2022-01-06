import ctypes as ct
import numpy as np 

class c_double_complex(ct.Structure):
	_fields_ = [("real", ct.c_double), ("imag", ct.c_double)]
	@property 
	def value(self):
		return self.real + 1j*self.imag 

def pointer_cmplx_array(cmplx_array):
	assert type(cmplx_array).__module__ == np.__name__, 'Input complex array must be numpy array'
	return cmplx_array.ctypes.data_as(ct.POINTER(c_double_complex))

def pointer_array(array):
	assert type(array).__module__ == np.__name__, 'Input array must be numpy array'
	return array.ctypes.data_as(ct.POINTER(ct.c_double))
