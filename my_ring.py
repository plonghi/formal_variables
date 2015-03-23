# import numpy as np
import sympy as sym
from sympy import poly, series # Using sympy v 0.7.6!
from sympy.abc import q, x, y, t

# FIND A WAY TO DISPENSE OF THIS GLOBAL VARIABLE EVENTUALLY
inverse_expansion_degree = 3

def dsz(charge_1, charge_2):
    return (charge_1[0] * charge_2[1] - charge_2[0] * charge_1[1])

class MyRing(object):
    """
    Elements of the ring are lists
    for example, r = x y + v w is represented by the list 
    [[x, y], [v, w]]
    Addition consists of list concatenation: if r' = a + b c then
    it is represented vt [[a], [b, c]] and r + r' is represented by
    [[x, y], [v, w], [a], [b, c]]
    Commutativity is not obvious at this level, but will be enforced 
    on the subclass representing the polynomials of formal variables.
    Multiplication is the distributed concatenation of the inner lists
    for example r * r' is represented by
    [[x, y, a], [x, y, b, c], [v, w, a], [v, w, b, c]]
    corresponding to
    x y a + x y b c + v w a + v w b c
    """
    
    def __init__(self, element):
        # make a check that the element is a listof depth 2
        self.element = element

    def __add__(self, other):
        el_1 = self.element
        el_2 = other.element
        new_el = el_1 + el_2 #concatenation
        return self.__class__(new_el) 
        # using self.__class__ is important for haivng sub-classes generate
        # other objects in the same sub-class.
        # if I used MyRing instead, the sum of two sub-classes would give
        # a MyRing type of object.

    def __mul__(self, other):
        el_1 = self.element
        el_2 = other.element
        new_el = [ \
                x + y for x in el_1 \
                for y in el_2
                ]
        return self.__class__(new_el)

    def __pow__(self, n):
        if type(n) is not int:
            print "Cannot elevate polynomial to a non-integer power."
            return None
        elif n == 0:
            return 1
        elif n < 0:
            if type(self) is not XVarRing:
                print "Cannot compute negative power of polynomial \
                that is not of type XVarRing (i.e. in formal variables)."
            else:
                c = self.laurent_filtered(0)  # the constant part
                # print "the constant element is: %s" % c
                if c == 0:
                    print "The polynomial %s vanishes at Y->0;\
                    \ncannot compute the inverse expansion" % self.__str__()
                    return None
                else:
                    ### e.g. if self is (3 + x)
                    norm_self = XVarRing([[ 1/c ]]) * self
                    # e.g. then norm_self is (1 + x/3)
                    # print norm_self
                    perturb_norm_self = norm_self + XVarRing([[-1]])
                    ### and then perturb_norm_self is + x/3
                    # print perturb_norm_self
                    inverse_norm_self = XVarRing([[1]])
                    for i in range(1, inverse_expansion_degree + 1):
                        inverse_norm_self = inverse_norm_self +\
                        XVarRing([[((-1) ** i)]])\
                        * (perturb_norm_self ** i)
                    ### then inverse_norm_self is 1 - x/3 + (x/3)^2 - (x/3)^3 
                    # print inverse_norm_self
                    inverse_self = XVarRing([[1/c]]) * inverse_norm_self
                    # print inverse_self
                    ### and inverse_self is the right answer for the inverse
                    negative_power_self = inverse_self ** abs(n)
                    return negative_power_self

        else:
            result = self.__class__([[1]])  # Multiplicative identity of ring
            for i in range(1, n+1):
                result = result * self
            return result



class FormalVariable(object):
    """
    The class of formal variables.
    """
    def __init__(self, charge):
        self.charge = charge

    def __str__(self):
        if self.charge == [0, 0]:
            return "1"
        else:
            return "Y_%s" % self.charge

    # def __mul__(self, other):
    #   ch_1 = self.charge
    #   ch_2 = other.charge
    #   new_ch = [ch_1[i] + ch_2[i] for i in range(len(ch_1))]
    #   return FormalVariable(new_ch)

    def laurent_form(self):
        return x**self.charge[0] * y**self.charge[1]



class XVarRing(MyRing):
    """
    This is the same as MyRing, but we use the convention that its entries 
    are of the form:
    [... [c, X_[q,p]] ...]
    where X_[q,p] is really an instance of the FormalVariable class.
    Actually the definition of this can be given in any form, and it will be 
    automatically rearranged in such a way
    """
    def __init__(self, element):
        super(XVarRing, self).__init__(element)
        self.simplify()

    def __str__(self):
        # self.simplify()
        el = self.element
        if el == [[]]:
            return "0"
        else:
            string = ""
            for j in range(len(self.element)):
                el = self.element[j]
                if el[0] == 1 and el[1].charge == [0, 0]:
                    new_piece = "1"
                elif el[1].charge == [0, 0] and not (el[0] == 1):
                    new_piece = " (%s) " % el[0]
                elif (not (el[1].charge == [0, 0])) and (el[0] == 1):
                    new_piece = el[1].__str__()
                else:
                    new_piece = (" (%s) " % el[0]) + el[1].__str__()
                string = string + new_piece
                if j < len(self.element) - 1: 
                    string = string + " + "
            return string

    def simplify(self):
        ### OLD
        # monomials = [twisted_product(x) for x in self.element]
        # self.element = map(de_laurentify, monomials)
        ### NEW
        monomials = [twisted_product(x) for x in self.element]
        # print "these are the monomials: %s" % monomials
        # Now I want to collect monomials with the same formal variable
        # and add together their coefficients
        ### OLD ###
        # self.element = resum(map(de_laurentify_m, monomials))
        ### NEW ###
        self.element = de_laurentify_p(sum(monomials))

    def represent(self):
        #### OLD
        # monomials = [twisted_product(x) for x in self.element]
        # return sum(monomials)
        # monomials = [twisted_product(x) for x in self.element]
        # self.simplify()
        ans = []
        for y in self.element:
            if len(y) == 0: 
                ans.append([])
            else: 
                ans.append([y[0], y[1].__str__()])
        return ans

    def laurentify_p(self):
        if self.element == [[]]:
            return 0
        else:
            laurent_poly = 0
            for e in self.element:
                summand = 1
                for f in e:
                    if type(f) is not FormalVariable:
                        summand *= f
                    else:
                        summand *= (x ** f.charge[0]) * (y ** f.charge[1])
                laurent_poly += summand
            return laurent_poly

    def t_expand(self, degree):
        laurent_poly = 0
        for e in self.element:
            summand = 1
            for f in e:
                if type(f) is not FormalVariable:
                    summand *= f
                else:
                    summand *= ((t * x) ** f.charge[0]) \
                            * ((t * y) ** f.charge[1])
            laurent_poly += summand
        return series(laurent_poly, t, 0, degree + 1).removeO()

    def laurent_filtered(self, degree):
        return self.t_expand(degree).subs(t, 1)

    def filtered(self, degree):
        return XVarRing(de_laurentify_p(self.laurent_filtered(degree)))





# def resum(l):
#     """
#     this is a function that takes
#     [... [c, X] , [d, Y] , [e, X] ...]
#     and 'resums it' into
#     [... [c + e, X] , [d, Y] ...]
#     it will also drop those terms such as 
#     [... [0, X] ...]
#     """
#     ### OLD
#     # var_form = [i[0] * i[1].laurent_form() for i in l]
#     # tot = sym.collect(sum(var_form), [x,y])
#     # term_list = tot.args
#     # print "\nthis is the resum function"
#     # print term_list
#     # return map(de_laurentify, term_list)

#     ### NEW 
#     # print "I want to resum this list: %s " % [[i[0], i[1].charge] for i in l]
#     # print "resumming %s " % l
#     tot = []
#     for i in l:
#         if len(i) == 0:
#             pass
#         else:
#             coeff = i[0]
#             var = i[1]
#             found = False
#             for k in range(len(tot)):
#                 # print "k = %s " % k
#                 if var.charge == tot[k][1].charge:
#                     # print "these are equal: %s = %s " % (var.charge , \
#                     # tot[k][1].charge)
#                     found = True
#                     position = k
#             if found == True:
#                 # print "will add %s to %s" % ( [coeff, var.charge], 
#                 # [tot[k][0], tot[k][1].charge] )
#                 tot[position][0] += coeff
#             elif found == False:
#                 tot.append(i)
#             # print "now tot is %s "% [ [s[0], s[1].charge] for s in tot]
#     ### Now will go through 'tot' and check if any cancellations occured, 
#     ### then will drop terms such as [...[0, X]...]
#     new_tot = []
#     for j in tot:
#         if j[0] == 0:
#             pass
#         else:
#             new_tot.append(j)
#     if len(new_tot) == 0:
#         new_tot = [[]]
#     # print "and returning %s " % tot
#     return new_tot


def twisted_product(l):
    """
    The twisted product between formal variables and other variables/numbers.
    For example, we associate X_[1,0] with the abelian variable x, 
    and X_[0,1] with y.
    Then we take
    [3, X_[a,b], X_[a',b'], 5]  
                    =>    15 * q^(<(a,b), (a',b')>) * x^(a+a') * y^(b+b')
    Note that the result is an expression in the abelina variables.
    There are other functions that will turn this back into formal variables,
    namely de_laurentify_m().
    """
    if len(l) == 0:
        return 0
    else:
        p = 1
        tot_charge = [0, 0]
        for i in l:
            if type(i) is FormalVariable:
                # print "this is a formal variable: %s" % i
                p *= i.laurent_form() * (q ** dsz(tot_charge, i.charge))
                tot_charge = [tot_charge[j] + i.charge[j] \
                                for j in range(len(i.charge))]
            else:
                p *= i
        return p


# def de_laurentify_m(m):
#     """
#     Takes an ABELIAN MONOMIAL in x,y variablea and translates it into an entry 
#     for an element of the Polynomial ring. For example
#         2 * q * x**a * y**b   =>   [2 * q, FormalVariable([a, b])]
#     """
#     # print "taking monomial %s "% m
#     if m == 0:
#         return [0, 0]
#     else:
#         d_x = sym.degree(m, x)
#         d_y = sym.degree(m, y)
#         c = m / x ** d_x / y ** d_y
#         # print type(m)
#         # print type(sym.poly(m, x, y))
#         # return (sym.poly(m, x, y)).subs(x,FormalVariable([1,0])).subs(\
#         #    y,FormalVariable([0,1]))
#         # return sym.poly(x).subs(x,FormalVariable([1,0])).subs(\
#         #    y,FormalVariable([0,1]))
#         # print "and returning %s "% [c, FormalVariable([d_x, d_y]).__str__()]
#         return [c, FormalVariable([d_x, d_y])]

def de_laurentify_p(p1):
    """
    Takes an ABELIAN POLYNOMIAL in x,y variablea and translates it into an 
    entry  for an element of the Polynomial ring. For example
        2 * q * x**a + 3 * y**b   
            =>   [[2 * q, FormalVariable([a, 0])], [3, FormalVariable([0, b])]]
    """
    # print "taking monomial %s "% m
    p = poly(p1, x, y)
    # print "type of %s is: %s" % (p1, type(p1))
    if p == 0:
        return [[]]
    else:
        ans = []
        coeffs = p.coeffs() ## the ordering is the same as in monoms()
        # print "coeffs: %s" % coeffs
        degrees = p.monoms() ## the ordering is lexicographic
        # print "monoms: %s" % degrees
        # d_x = sym.degree(m, x)
        # d_y = sym.degree(m, y)
        # c = m / x ** d_x / y ** d_y
        # print type(m)
        # print type(sym.poly(m, x, y))
        # return (sym.poly(m, x, y)).subs(x,FormalVariable([1,0])).subs(\
        #    y,FormalVariable([0,1]))
        # return sym.poly(x).subs(x,FormalVariable([1,0])).subs(\
        #    y,FormalVariable([0,1]))
        # print "and returning %s "% [c, FormalVariable([d_x, d_y]).__str__()]
        for i in range(len(coeffs)):
            ans.append([coeffs[i], FormalVariable(list(degrees[i]))])
        return ans



# e1 = MyRing([[1,2]])
# e2 = MyRing([[3, 4, 5]])
# e3 = MyRing([[x, y]])
# print e1.element
# print e2.element
# print ( e1 + e2 ).element
# print ( e1 * e2 ).element
# print ( (e1 + e2) * e3 ).element
# print ( e1 * e3 + e2 * e3 ).element
# print (e1 + (e2 + e3)).element
# print ((e1 + e2) + e3).element
# print (e1 * (e2 * e3)).element
# print ((e1 * e2) * e3).element


init_message = \
"""
To define a formal variable with charge [p, q]
use the command:
    x1 = FormalVariable([p, q])
To define an element of the ring of formal variables
use the command:
    P = XVarRing([ ..., [x1, x2, 3, q ], ...  ])
for a polynomial containing
    ... + 3 q x1 x2 + ...
Note that the formal variables don't commute in general,
while integers and the quantum parameter q are abelian.
"""

print init_message

#-------------------------------#
#            EXAMPLE            #
#-------------------------------#

# print "\nan abstract ring element"
# e4 = MyRing([[2,x], [y], [3, x, y, x]])
# print e4.element
# print "\nits square"
# print type(e4 **2)
# print (e4 ** 2).element

print "\nlet's define some formal variables"
x1 = FormalVariable([1,0])
x2 = FormalVariable([0,1])
print "x1 = %s" % x1
print "x2 = %s" % x2

print "\nlet's define some polynomials"
e5 = XVarRing([[x1, x2], [x2]])
print "e5 = %s" % e5
e6 = XVarRing([[x2, x1], [x2]])
print "e6 = %s" % e6

print "\nlet's multiply the polynomials"
print "e5 * e6 = %s" % (e5 * e6)

print "\nlet's add the polynomials"
print "e5 + e6 = %s" % (e5 + e6)

print "\nlet's take the square power of e5"
print "e5^2 = %s" % e5 ** 2

print "\nlet's take the inverse square power of e5"
print "e5^{-2} = %s" % e5 ** (-2)

print "\nlet's try with the inverse square power of e7"
# e7 = XVarRing([[1], [2, q, x1], [3, 1 / q**2, x2, x1]])
e7 = XVarRing([[1], [x1], [x2]])
# e7 = XVarRing([[1], [q, x1]])
# e7 = XVarRing([[3], [x1]])
print "e7 = %s" % e7 
print "e7^{-2} = %s" % e7 ** (-2)
print "\na 'filtered' presentation of this expression"
print "e7^{-2} = %s" % (e7 ** (-2)).filtered(inverse_expansion_degree)

# print "\nThe 'laurent abelianized' expression of e5:"
# print e5.laurentify_p()

# k = 1
# print "\nThe 'laurent abelianized' t-expanded expression of e5 up to \
# degree t^%s: " % k
# print e5.t_expand(k)

# print "\nThe filtered truncation of e5 up to degree %s:" % k
# print e5.laurent_filtered(k)

# print "\nThe 'laurent abelianized' expression of e5 + e5 * e6:"
# f1 = (e5 + e5 * e6).laurentify_p()
# print f1
# print "\nThe 'de-laurentified' expression of e5 + e5 * e6:"
# print XVarRing(de_laurentify_p(f1))

# print "\n\nSome special cases:"

# print "\nDefining \
# \ne8 = XVarRing([[1], [-1]])\
# \nthe sum of additive inverses"
# e8 = XVarRing([[1], [-1]])
# print "after automatic simplification, it has 'elements' %s " % e8.element
# print "e8 = %s" % e8
# print "the laurent-abelianized-then de_laurentified expression: %s " % \
# XVarRing(de_laurentify_p(e8.laurentify_p()))

# print "\nDefining\
# \ne9 = XVarRing([[0]])\
# \nthe additive identity"
# e9 = XVarRing([[0]])
# print "after automatic simplification, it has 'elements' %s " % e9.element
# print "e9 = %s" % e9
# print "the laurent-abelianized-then de_laurentified expression: %s " % \
# XVarRing(de_laurentify_p(e9.laurentify_p()))

# print "\nDefining\
# \ne10 = XVarRing([[1]])\
# \ntha multiplicative identity"
# e10 = XVarRing([[1]])
# print "after automatic simplification, it has 'elements' %s " % e10.element
# print "e10 = %s" % e10
# print "the laurent-abelianized-then de_laurentified expression: %s " % \
# XVarRing(de_laurentify_p(e10.laurentify_p()))

# print "\nDefining \
# \ne11 = XVarRing([[FormalVariable([1, 2]), \
# FormalVariable([-1, -2])]])\
# \nthe product of multiplicative inverses"
# e11 = XVarRing([[FormalVariable([1, 2]), FormalVariable([-1, -2])]])
# print "after automatic simplification, it has 'elements' %s " % e11.element
# print "e11 = %s" % e11
# print "the laurent-abelianized-then de_laurentified expression: %s " % \
# XVarRing(de_laurentify_p(e11.laurentify_p()))

#-------------------------------#
#         END EXAMPLE           #
#-------------------------------#

# print twisted_product([ x1, x1 ])
# print twisted_product([ x1, x2 ])
# print de_laurentify( twisted_product([ x2, x1 ]) )

# r = twisted_product([ x1, x2 ])
# print de_laurentify( x )
# print de_laurentify( x )[1]

# e7 = e5 * e6
# print e7.element
# print e7.represent()