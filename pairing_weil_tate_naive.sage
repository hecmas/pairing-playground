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
def line(P, Q, D):
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
        f = f * f * line(R,R,D) / line(-2*R,2*R,D)
        R = 2 * R
        if r & 2^i:
            f = f * line(R,P,D) / line(-(R+P),R+P,D)
            R = R + P

    f = f * f * line(R,R,D)

    return f

def untwist(eE,Q,d):
    Fpk = eE.base_field()

    if Q.is_zero():
        return Q

    assert d in [2,3,4,6]
    # The embedding degree k should be divisible by d
    assert k % d == 0

    if d == 6 or d == 3:
        # E should be of the form y² = x³ + b
        assert eE.a4() == 0
    elif d == 4:
        # E should be of the form y² = x³ + ax
        assert eE.a6() == 0

    _x,_y = Q.xy()
    # Isomorphism from F[i]/<i²+1> to F[j]/<j² - 2·j + 2>
    xcoeffs = [_x.polynomial().list()[0] - _x.polynomial().list()[1],_x.polynomial().list()[1]]
    ycoeffs = [_y.polynomial().list()[0] - _y.polynomial().list()[1],_y.polynomial().list()[1]]

    # Isomorphism into the subfield of F[w]/<w¹² - 2·w⁶ + 2>,
    # where w⁶ = j
    nx = Fpk(xcoeffs[0] + w^6*xcoeffs[1])
    ny = Fpk(ycoeffs[0] + w^6*ycoeffs[1])

    return eE(nx / w^2, ny / w^3)

def Weil(P,Q):
    if P.is_zero() or Q.is_zero():
        return F.one()
    
    DP = Divisor([1,-1],[2*P,P])
    DQ = Divisor([1,-1],[2*Q,Q])
    frP = Miller_Loop(P,DQ)
    frQ = Miller_Loop(Q,DP)
    lPP = line(P,P,DQ)
    v2P = line(-2*P,2*P,DQ)
    lQQ = line(Q,Q,DP)
    v2Q = line(-2*Q,2*Q,DP)
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
# q = 19
# F = GF(q)
# E = EllipticCurve(F, [14,3])
# n = E.order()
# # the torsion group parameter r is typically chosen 
# # as the largest prime factor of the order of the curve
# r = n.factor()[-1][0]
# k = embedding_degree(q,r)

# P = E(17,9)
# assert r*P == E(0) # P is in the r-torsion subgroup

# # To define Q, we need to move to the extension field F_{q**k}
# K.<i> = GF(q^k, modulus=x^2+1)
# PRK.<x,y> = PolynomialRing(K)
# eE = EllipticCurve(K, [14,3])

# Q = eE(16, 16*i)
# assert r*Q == eE(0) # Q is in the r-torsion subgroup

# assert Tate(P,Q) == 15*i + 2

# # For sure, the pairing is non-degenerate
# assert P.additive_order() == Q.additive_order() == r
# assert Tate(P,Q) != F.one()

# # Let's check the bilinearity of the pairing
# assert Tate(2*P,2*Q) == Tate(P,2*Q)^2 == Tate(2*P,Q)^2 == Tate(P,Q)^4 == Tate(2*P,2*Q)

# # Check the trivial evaluations are satisfied
# assert Tate(E(0),Q) == Tate(P,E(0)) == F.one()

# # Since P and Q are generators, we should have that Tate(P,D) is a primitive r-th root of unity
# # i.e. a generator of set of roots of unity of order r
# assert multiplicative_order(Tate(P,Q)) == r

# # Test 3: Tate Pairing over the BLS12-381 curve
# # https://hackmd.io/@benjaminion/bls12-381
# bls_x = -15132376222941642752
# q = 1/3*(bls_x-1)^2*(bls_x^4 - bls_x^2 + 1) + bls_x
# assert len(bin(q)[2:]) == 381
# r = bls_x^4 - bls_x^2 + 1
# k = embedding_degree(q,r)
# assert 12 == k
# F = GF(q)
# E = EllipticCurve(F, [0,4])

# P = E(3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507,
#       1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569)
# # Subgrup check: P is in G1 = E(Fp)[r]
# assert r*P == E(0)

# # To define Q, we need to move to the extension field Fp2
# Fp2.<i> = GF(q^2, modulus=x^2+1)
# tE = EllipticCurve(Fp2, [0,4*(i+1)])

# Q = tE(Fp2(3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758*i 
#         + 352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160),
#         Fp2(927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582*i
#         + 1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905))

# # Subgrup check: Q is in G2' = E'(Fp2)[r]
# assert r*Q == tE(0) # Q is in the r-torsion subgroup

# Fp12.<w> = GF(q^12, modulus=x^12 - 2*x^6 + 2)
# PRK.<x,y> = PolynomialRing(Fp12)
# eE = E.base_extend(Fp12)

# # Less efficient subgrup check: Q is in G2 = E(Fp12)[r] ∩ Ker(Frob - [q])
# tQ = untwist(eE,Q,6)
# assert r*tQ == eE(0)
# FrobtQ = eE(tQ[0]^q,tQ[1]^q)
# QQ = int(q)*tQ
# assert FrobtQ == QQ

# # For sure, the pairing is non-degenerate
# # assert P.additive_order() == tQ.additive_order() == r
# assert Tate(P,untwist(eE,Q,6)) != F.one()

# # Let's check the bilinearity of the pairing
# assert Tate(2*P,untwist(eE,2*Q,6)) == Tate(P,untwist(eE,2*Q,6))^2 == Tate(2*P,untwist(eE,Q,6))^2 == Tate(P,untwist(eE,Q,6))^4 == Tate(2*P,untwist(eE,2*Q,6))

# # Check the trivial evaluations are satisfied
# assert Tate(E(0),untwist(eE,Q,6)) == Tate(P,untwist(eE,tE(0),6)) == F.one()

# # Since P and Q are generators, we should have that Tate(P,D) is a primitive r-th root of unity
# # i.e. a generator of set of roots of unity of order r
# # assert multiplicative_order(Tate(P,untwist(eE,Q,6))) == r