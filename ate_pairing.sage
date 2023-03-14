import os
os.system('sage --preparse tools.sage')
os.system('mv tools.sage.py tools.py')
from tools import log2, embedding_degree, line

# This is the ate pairing obtained with loop shortening optimizations
# Note: It only works for even embedding degrees k
# r is assumed to be in binary form as r[i]
# For sure, Q is assumed to be from the (only) subgroup of E[r] over F_q
# P is assumed to be from the trace zero subgroup of E[r] over F_{q**k}
# i.e., a (reversed) Type 3 pairing is assumed
def Miller_Loop(Q,P):
    if P.is_zero() == True or Q.is_zero() == True:
        return F.one()

    R = Q
    f = F.one()
    for i in range(log2(T)-2,-1,-1):
        f = f * f * line(R,R,P,E)
        R = 2 * R
        if T & 2^i:
            f = f * line(R,Q,P,E)
            R = R + Q

    return f

def ate(Q,P):
    return Miller_Loop(Q,P)^((q^k-1)/r)

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

t = q+1-n
T = abs(t-1)
assert ate(Q,P) == 21*u^3 + 37*u^2 + 25*u + 25

# For sure, the pairing is non-degenerate
assert P.additive_order() == Q.additive_order() == r
assert ate(Q,P) != F.one()

# Let's check the bilinearity of the pairing
assert ate(2*Q, 12*P) == ate(Q,12*P)^2 == ate(2*Q,P)^12 == ate(Q,P)^24 == ate(12*Q,2*P)

# Check the trivial evaluations are satisfied
assert ate(E(0),P) == ate(Q,E(0)) == F.one()

# Since P and Q are generators, we should have that ate(Q,P) is a primitive r-th root of unity
# i.e. a generator of set of roots of unity of order r
assert multiplicative_order(ate(Q,P)) == r