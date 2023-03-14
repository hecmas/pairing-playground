import os
os.system('sage --preparse tools.sage')
os.system('mv tools.sage.py tools.py')
from tools import log2, embedding_degree, line

# This is the BKLS-GHS version of Miller's algorithm for computing the Weil and Tate pairings
# Note: It only works for even embedding degrees k
# r is assumed to be in binary form as r[i]
# For sure, P is assumed to be from the (only) subgroup of E[r] over F_q 
# Q is assumed to be from the trace zero subgroup of E[r] over F_{q**k}
# i.e., a Type 3 pairing is assumed
def Miller_Loop(P,Q):
    if P.is_zero() == True or Q.is_zero() == True:
        return F.one()

    R = P
    f = F.one()
    # Only loop until the second to last bit of r
    for i in range(log2(r)-2,0,-1):
        f = f * f * line(R,R,Q,E)
        R = 2 * R
        if r & 2^i:
            f = f * line(R,P,Q,E)
            R = R + P

    f = f * f * line(R,R,Q,E)

    return f

def Tate(P,Q):
    return Miller_Loop(P,Q)^((q^k-1)/r)

# Weil is not working, fix it!
def Weil(P,Q):
    a = Miller_Loop(P,Q)
    b = Miller_Loop(Q,P)
    print(a,b)
    return (-1)^r * a / b

q = 47
F = GF(q)
E = EllipticCurve(F, [21,15])
n = E.order()
# the torsion group parameter r is typically chosen 
# as the largest prime factor of the order of the curve
r = n.factor()[-1][0]
k = embedding_degree(q,r)
assert k % 2 == 0

P = E(45,23)
assert r*P == E(0) # P is in the r-torsion subgroup

# To define Q, we need to move to the extension field F_{q**k}
K.<u> = GF(q^k, modulus=x^4-4*x^2+5)
eE = EllipticCurve(K, [21,15])

Q = eE(31*u^2 + 29, 35*u^3 + 11*u)
assert 17*Q == eE(0) # Q is in the r-torsion subgroup
#TODO: Check that Q is in the trace zero subgroup

assert Tate(P,Q) == 33*u^3 + 43*u^2 + 45*u + 39

# For sure, the pairing is non-degenerate
assert P.additive_order() == Q.additive_order() == r
assert Tate(P,Q) != F.one()

# Let's check the bilinearity of the pairing
assert Tate(2*P, 12*Q) == Tate(P,12*Q)^2 == Tate(2*P,Q)^12 == Tate(P,Q)^24 == Tate(12*P,2*Q)
# assert Weil(2*P, 12*Q) == Weil(P,12*Q)^2 == Weil(2*P,Q)^12 == Weil(P,Q)^24 == Weil(12*P,2*Q)

# Check the trivial evaluations are satisfied
assert Tate(E(0),Q) == Tate(P,E(0)) == F.one()

# Since P and Q are generators, we should have that Tate(P,Q) is a primitive r-th root of unity
# i.e. a generator of set of roots of unity of order r
assert multiplicative_order(Tate(P,Q)) == r