def embedding_degree(E):
    # n = E.order()
    # r = n.factor()[-1][0]
    
    # we know that |E| = q+1-t, where t is the trace of Frobenius
    t = E.trace_of_frobenius()
    q = n - 1 + t

    k = 1
    while (q**k - 1) % r != 0:
        k += 1

    return k

# Find line y = mx + c passing through two points P and Q
# or vertical line y = x0 if Q = -P
# and evaluate it at a point T
def l(P, Q, T):
    assert P.is_zero() != True and Q.is_zero() != True and T.is_zero() != True

    # First case: P and Q are distinct and not on the same vertical line
    if P.xy()[0] != Q.xy()[0]:
        m = (Q.xy()[1] - P.xy()[1]) / (Q.xy()[0] - P.xy()[0])
        c = P.xy()[1] - m * P.xy()[0]
        return T.xy()[1] - m * T.xy()[0] - c
    # Second case: P and Q are the same point
    elif P.xy()[1] == Q.xy()[1]:
        m = (3 * P.xy()[0] * P.xy()[0] + E.a4()) / (2 * P.xy()[1])
        c = P.xy()[1] - m * P.xy()[0]
        return T.xy()[1] - m * T.xy()[0] - c
    # Third case: P and Q are distinct and on the same vertical line
    # The line is y = P.xy()[0]
    else:
        return T.xy()[1] - P.xy()[0]

# This is the BKLS-GHS version of Miller's algorithm for computing the ate pairings
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
    for i in range(len(Tbin)-2,-1,-1):
        f = f * f * l(R,R,P)
        R = 2 * R
        if Tbin[i] == 1:
            f = f * l(R,Q,P)
            R = R + Q

    # Here, k is the embedding degree of E
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
k = embedding_degree(E)
assert k % 2 == 0

P = E(45,23)
assert r*P == E(0) # P is in the r-torsion subgroup

# To define Q, we need to move to the extension field F_{q**k}
K.<u> = GF(q^k, modulus=x^4-4*x^2+5)
eE = EllipticCurve(K, [21,15])

Q = eE(31*u^2 + 29, 35*u^3 + 11*u)
assert 17*Q == eE(0) # Q is in the r-torsion subgroup

t = E.trace_of_frobenius()
T = abs(t-1)
Tbin = T.digits(2)
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