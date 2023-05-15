import os
os.system('sage --preparse divisors.sage')
os.system('mv divisors.sage.py divisors.py')
from divisors import evaluate_function_on_divisor,Divisor

# Test 1: Divisor playground
q = 47
F = GF(q)
E = EllipticCurve(F, [21,15])
D1 = Divisor(E, [2, -3], [E(1,15), E(1,32)])
D2 = Divisor(E, [3, 1, -1], [E(1,32), E(2,21), E(6,13)])
D3 = D1 + D2
assert D1.degree == -1
assert D2.degree == 3
assert D3 == Divisor(E, [2, 1, -1], [E(1,15), E(2,21), E(6,13)])
assert D3.degree == D1.degree + D2.degree == 2
assert D1.support == [E(1,15), E(1,32)]
assert D2.support == [E(1,32), E(2,21), E(6,13)]
assert D3.support == [E(1,15), E(2,21), E(6,13)]

# The zero divisor
D4 = Divisor.zero(E)
D5 = Divisor(E,[],[])
assert (D1 - D1).is_zero() == True and D4 == D5 and (D1 - D1) == D4

D = Divisor(E, [1,1,-2], [E(1,15),E(1,32),E(0)])
assert D.is_principal() == True


# Test 2: Divisor equivalence
q = 61
F = GF(q)
E = EllipticCurve(F, [8,1])
P = E(57,24)
Q = E(25,37)
R = E(17,32)
S = E(42,35)
D1 = Divisor(E,[1,1,1], [P,Q,R])
D2 = Divisor(E,[4,-1], [E(0),S])
assert D1.equivalent(D2) == True


# Test 3: Function evaluation on divisor
q = 163
F = GF(q)
E = EllipticCurve(F, [-1,-2])
P = E(43,154)
Q = E(46,38)
R = E(12,35)
S = E(5,66)

PR.<x,y> = PolynomialRing(F)
lPQ = y + 93*x + 85
lPP = y + 127*x + 90
lQQ = y + 13*x + 16
assert lPQ(P.xy()) == lPQ(Q.xy()) == 0
assert lPP(P.xy()) == lPP((-2*P).xy()) == 0
assert lQQ(Q.xy()) == lQQ((-2*Q).xy()) == 0
D1 = Divisor(E,[2,1],[R,S])
D2 = Divisor(E,[3,-3],[R,S])
D3 = Divisor(E,[1,1,-2],[R,S,E(0)])
assert evaluate_function_on_divisor(lPQ,D1) == 122
assert evaluate_function_on_divisor(lPP,D2) == 53

# Test 4: Weil Reciprocity
q = 503
F = GF(q)
E = EllipticCurve(F, [0,1])
P1,P2,P3,P4 = E(433,98),E(232,113),E(432,27),E(127,258)
Q1,Q2,Q3,Q4 = E(413,369),E(339,199),E(147,443),E(124,42)

# TODO: Implement f --> (f)

PR.<x,y> = PolynomialRing(F)
f = (20*y + 9*x + 179)/(199*y + 187*x + 359)
g = y + 251*x^2 + 129*x + 201
# Here, Df = (f) and Dg = (g)
Df = Divisor(E,[2,1,-1,-2],[P1,P2,P3,P4])
Dg = Divisor(E,[1,1,1,1,-4],[Q1,Q2,Q3,Q4,E(0)])
# If (f) and (g) have disjoint supports, then Weil reciprocity ensures that f((g)) = g((f))
assert evaluate_function_on_divisor(f,Dg) == evaluate_function_on_divisor(g,Df) == 321
