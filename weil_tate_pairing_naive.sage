import os
os.system('sage --preparse divisors.sage')
os.system('mv divisors.sage.py divisors.py')
os.system('sage --preparse tools.sage')
os.system('mv tools.sage.py tools.py')
from divisors import evaluate_function_on_divisor,Divisor
from tools import log2, embedding_degree

# Find line y = mx + c passing through two points P and Q
# or vertical line y = x0 if Q = -P
# and evaluate it at a divisor D
def l(P, Q, D):
    assert P.is_zero() != True and Q.is_zero() != True

    # First case: P and Q are distinct and not on the same vertical line
    if P.xy()[0] != Q.xy()[0]:
        m = (Q.xy()[1] - P.xy()[1]) / (Q.xy()[0] - P.xy()[0])
        c = P.xy()[1] - m * P.xy()[0]

        f = y - m * x - c
        return evaluate_function_on_divisor(f,D)
    # Second case: P and Q are the same point
    elif P.xy()[1] == Q.xy()[1]:
        m = (3 * P.xy()[0] * P.xy()[0] + E.a4()) / (2 * P.xy()[1])
        c = P.xy()[1] - m * P.xy()[0]
        f = y - m * x - c
        return evaluate_function_on_divisor(f,D)
    # Third case: P and Q are distinct and on the same vertical line
    # The line is x = P.xy()[0]
    else:
        f = x - P.xy()[0]
        return evaluate_function_on_divisor(f,D)

# This is the naive version of Miller's algorithm for computing the Weil and Tate pairings
# For sure, P is assumed to be from the (only) subgroup of E[r] over F_q 
# Q is assumed to be from the trace zero subgroup of E[r] over F_{q**k}
# i.e., a Type 3 pairing is assumed
def Miller_Loop(P,D):
    R = P
    f = F.one()
    # Only loop until the second to last bit of r
    for i in range(log2(r)-2,0,-1):
        f = f * f * l(R,R,D) / l(-2*R,2*R,D)
        R = 2 * R
        if r & 2^i:
            f = f * l(R,P,D) / l(-(R+P),R+P,D)
            R = R + P

    f = f * f * l(R,R,D)

    return f

def Weil(P,Q):
    if P.is_zero() or Q.is_zero():
        return F.one()
    
    DP = Divisor([1,-1],[2*P,P])
    DQ = Divisor([1,-1],[2*Q,Q])
    frP = Miller_Loop(P,DQ)
    frQ = Miller_Loop(Q,DP)
    lPP = l(P,P,DQ)
    v2P = l(-2*P,2*P,DQ)
    lQQ = l(Q,Q,DP)
    v2Q = l(-2*Q,2*Q,DP)
    A = frP / (lPP / v2P)^3
    B = frQ / (lQQ / v2Q)^3
    return A/B

def Tate(P,Q):
    if P.is_zero() or Q.is_zero():
        return F.one()
    
    DQ = Divisor([1,-1],[2*Q,Q])

    # Here, k is the embedding degree of E
    return Miller_Loop(P,DQ)^((q^k-1)/r)

# Test 1: Weil Pairing
q = 23
F = GF(q)
E = EllipticCurve(F, [-1,0])
n = E.order()
# the torsion group parameter r is typically chosen 
# as the largest prime factor of the order of the curve
r = n.factor()[-1][0]
k = embedding_degree(q,r)

P = E(2,11)
assert r*P == E(0) # P is in the r-torsion subgroup

# To define Q, we need to move to the extension field F_{q**k}
K.<i> = GF(q^k, modulus=x^2+1)
PRK.<x,y> = PolynomialRing(K)
eE = EllipticCurve(K, [-1,0])

Q = eE(21, 12*i)
assert r*Q == eE(0) # Q is in the r-torsion subgroup

assert Weil(P,Q) == 15*i + 11

# For sure, the pairing is non-degenerate
assert P.additive_order() == Q.additive_order() == r
assert Weil(P,Q) != F.one()

# Let's check the bilinearity of the pairing
assert Weil(2*P,2*Q) == Weil(P,2*Q)^2 == Weil(2*P,Q)^2 == Weil(P,Q)^4 == Weil(2*P,2*Q)

# Check the trivial evaluations are satisfied
assert Weil(E(0),Q) == Weil(P,E(0)) == F.one()

# Since P and Q are generators, we should have that Weil(P,D) is a primitive r-th root of unity
# i.e. a generator of set of roots of unity of order r
assert multiplicative_order(Weil(P,Q)) == r


# Test 2: Tate Pairing
q = 19
F = GF(q)
E = EllipticCurve(F, [14,3])
n = E.order()
# the torsion group parameter r is typically chosen 
# as the largest prime factor of the order of the curve
r = n.factor()[-1][0]
k = embedding_degree(q,r)

P = E(17,9)
assert r*P == E(0) # P is in the r-torsion subgroup

# To define Q, we need to move to the extension field F_{q**k}
K.<i> = GF(q^k, modulus=x^2+1)
PRK.<x,y> = PolynomialRing(K)
eE = EllipticCurve(K, [14,3])

Q = eE(16, 16*i)
assert r*Q == eE(0) # Q is in the r-torsion subgroup

assert Tate(P,Q) == 15*i + 2

# For sure, the pairing is non-degenerate
assert P.additive_order() == Q.additive_order() == r
assert Tate(P,Q) != F.one()

# Let's check the bilinearity of the pairing
assert Tate(2*P,2*Q) == Tate(P,2*Q)^2 == Tate(2*P,Q)^2 == Tate(P,Q)^4 == Tate(2*P,2*Q)

# Check the trivial evaluations are satisfied
assert Tate(E(0),Q) == Tate(P,E(0)) == F.one()

# Since P and Q are generators, we should have that Tate(P,D) is a primitive r-th root of unity
# i.e. a generator of set of roots of unity of order r
assert multiplicative_order(Tate(P,Q)) == r