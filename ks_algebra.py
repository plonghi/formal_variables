# import numpy as np
import math
import my_ring
from my_ring import XVarRing as P, FormalVariable as Y
from misc import k_delta
import sympy as sym
from sympy import poly # Using sympy v 0.7.6!
from sympy.abc import x, y, q

lattice_rank = 2
### to be generalized

basis_vectors = [\
                [k_delta(i,j)  for j in range(lattice_rank)] \
                for i in range(lattice_rank)]
# print basis_vectors

basis_variables = [P([[Y(v)]]) for v in basis_vectors]
### Note: this is a list of RING ELEMENTS, each is a POLY-(really MONO-)NOMIAL 
print "\nThe basis of variables is:"
for p in basis_variables:
    print p


def dilog_factor(s, sign, ring_monomial):
    coeff = ring_monomial.element[0][0]
    var = ring_monomial.element[0][1]
    return P([[1], [q**((-1) * (sign) * (2 * s -1)) * coeff, var]])

# print dilog_factor(1, 1, basis_variables[0])


def dilog(n, ring_m):
    """
    Finite type quantum dilog. The ring_m should be a monomial of the ring,
    hence of the form [[2, Y_[a,b]]]
    """
    if not (type(ring_m) is my_ring.XVarRing and type(n) is int):
        print "Wrong variable type in Dilogarithm."
        return None
    else:
        pass
    
    if n == 0:
        return 1
    elif n > 0:
        sign = 1
    else:
        sign = -1

    result = P([[1]])  # The identity in the polynomial ring
    for i in range(1, abs(n)+1):
        result = result * dilog_factor(i, sign, ring_m)

    return result


#-------------------------------#
#            TO DO              #
#-------------------------------#

def KS(ring_P, PSC):
    if type(ring_P) is my_ring.XVarRing:
        print "good."
    else:
        print "%s is not a valid element of the polynomial ring" % ring_P




print basis_variables[0] \
        * dilog(1, P([[ (-1 * q)**(2) ]]) * basis_variables[1]) \
        * dilog(1, P([[ (-1 * q)**(0) ]]) * basis_variables[1]) \
        * dilog(1, P([[ (-1 * q)**(-2) ]]) * basis_variables[1]) \