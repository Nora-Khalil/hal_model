import cantera as ct
import numpy as np
import pandas as pd

print("Running Cantera Version: " + str(ct.__version__))

To = 298
Po = 1e5 # ct.one_atm

gas = ct.Solution('./cantera/chem.cti')


bdc = gas.binary_diff_coeffs
bdcpoly = gas.get_binary_diff_coeffs_polynomial(1,3)

print('binary_diff_coeffs are:')
print(bdc) 
print('binary_diff_coeffs_polynomials are:')
print(bdcpoly)
