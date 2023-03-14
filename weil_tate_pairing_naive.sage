import os
os.system('sage --preparse divisors.sage')
os.system('mv divisors.sage.py divisors.py')
from divisors import evaluate_function_on_divisor,Divisor

def embedding_degree(q,r):    
    k = 1
    while (q**k - 1) % r != 0:
        k += 1

    return k

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

# This is the BKLS-GHS version of Miller's algorithm for computing the Weil and Tate pairings
# Note: It only works for even embedding degrees k
# r is assumed to be in binary form as r[i]
# For sure, P is assumed to be from the (only) subgroup of E[r] over F_q 
# Q is assumed to be from the trace zero subgroup of E[r] over F_{q**k}
# i.e., a Type 3 pairing is assumed
def Miller_Loop(P,D):
    if P.is_zero() == True or D == Divisor.zero():
        return F.one()

    R = P
    f = F.one()
    # Only loop until the second to last bit of r
    for i in range(len(rbin)-2,0,-1):
        f = f * f * l(R,R,D) / l(-2*R,2*R,D)
        R = 2 * R
        if rbin[i] == 1:
            f = f * l(R,P,D) / l(-(R+P),R+P,D)
            R = R + P

    f = f * f * l(R,R,D)

    return f

def Weil(P,Q):
    DP = Divisor([1,-1],[2*P,P])
    DQ = Divisor([1,-1],[2*Q,Q])
    a = Miller_Loop(P,DQ)
    b = Miller_Loop(Q,DP)
    c = (l(Q,Q,DP) / l(-2*Q,2*Q,DP))^r
    print("a",a)
    print("b",b)
    print("c",c)
    return a/(b/c)

def Tate(P,Q):
    DQ = Divisor([1,-1],[2*Q,Q])

    # Here, k is the embedding degree of E
    return Miller_Loop(P,DQ)^((q^k-1)/r)

q = 47
F = GF(q)
PR.<x,y> = PolynomialRing(F)
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

rbin = r.digits(2)
# assert Tate(P,Q) == 33*u^3 + 43*u^2 + 45*u + 39
print(Weil(P,Q))
# assert Weil(P,D) == 22*u^3 + 12*u^2 + 32*u + 13

# For sure, the pairing is non-degenerate
assert P.additive_order() == Q.additive_order() == r
# assert Weil(P,Q) != F.one()

# # Let's check the bilinearity of the pairing
assert Tate(2*P, 12*Q) == Tate(P,12*Q)^2 == Tate(2*P,Q)^12 == Tate(P,Q)^24 == Tate(12*P,2*Q)
# assert Weil(2*P, 12*Q) == Weil(P,12*Q)^2 == Weil(2*P,Q)^12 == Weil(P,Q)^24 == Weil(12*P,2*Q)

# # Check the trivial evaluations are satisfied
# assert Weil(E(0),Q) == Weil(P,E(0)) == F.one()

# # Since P and Q are generators, we should have that Weil(P,D) is a primitive r-th root of unity
# # i.e. a generator of set of roots of unity of order r
# assert multiplicative_order(Weil(P,Q)) == r