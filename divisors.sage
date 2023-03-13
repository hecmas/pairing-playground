
# TODO: Add the zero divisor
class Divisor():
    def __init__(self, multiplicites, points):
        assert len(multiplicites) == len(points)
        for m in multiplicites:
            assert m in ZZ
        # This is not necessary, assuming that every point P in points is defined as E(P)
        # for P in points:
        #     assert P in E

        self.multiplicites = multiplicites
        self.points = points
        self.degree = sum(multiplicites)
        self.support = points

        self.sum = {}
        for i in range(len(multiplicites)):
            self.sum[points[i]] = multiplicites[i]

    def __repr__(self):
        result = ""
        for i in range(len(self.multiplicites)-1):
            result += str(self.multiplicites[i]) + " " + str(self.points[i]) + " + "

        result += str(self.multiplicites[-1]) + " " + str(self.points[-1])
        
        return result

    def __add__(self, other):
        multiplicites = []
        points = []
        for key in self.sum.keys():
            if key in other.sum.keys():
                keysum = self.sum[key] + other.sum[key]
                if keysum != 0:
                    multiplicites.append(keysum)
                    points.append(key)
            else:
                multiplicites.append(self.sum[key])
                points.append(key)

        for key in other.sum.keys():
            if key not in self.sum.keys():
                multiplicites.append(other.sum[key])
                points.append(key)

        return Divisor(multiplicites, points)
    
    def __neg__(self):
        for key in self.sum.keys():
            self.sum[key] = -self.sum[key]
        return self

    def __sub__(self, other):
        return self + (-other)
    
    def __eq__(self, other):
        return self.sum == other.sum
    
    def __ne__(self, other):
        return not self == other
    
    def is_principal(self):
        addition = E(0)
        for i in range(len(self.multiplicites)):
            addition += self.multiplicites[i] * E(self.points[i])

        if addition == E(0) and self.degree == 0:
            return True
        else:
            return False

def evaluate_function_on_divisor(f,D):
    # TODO: Add a check for disjoint supports
    result = 1
    for i in range(len(D.multiplicites)):
        result *= f(D.points[i].xy()) ^ D.multiplicites[i]
    return result

# Test 1
# q = 47
# F = GF(q)
# E = EllipticCurve(F, [21,15])
# D1 = Divisor([2, -3], [E(1,15), E(1,32)])
# D2 = Divisor([3, 1, -1], [E(1,32), E(2,21), E(6,13)])
# D3 = D1 + D2
# assert D1.degree == -1
# assert D2.degree == 3
# assert D3 == Divisor([2, 1, -1], [E(1,15), E(2,21), E(6,13)])
# assert D3.degree == D1.degree + D2.degree == 2
# assert D1.support == [E(1,15), E(1,32)]
# assert D2.support == [E(1,32), E(2,21), E(6,13)]
# assert D3.support == [E(1,15), E(2,21), E(6,13)]

# D = Divisor([1,1,-2],[E(1,15),E(1,32),E(0)])
# assert D.is_principal() == True

# Test 2
# q = 163
# F = GF(q)
# E = EllipticCurve(F, [-1,-2])
# P = E(43,154)
# Q = E(46,38)
# R = E(12,35)
# S = E(5,66)

# PR.<x,y> = PolynomialRing(F)
# lPQ = y + 93*x + 85
# lPP = y + 127*x + 90
# lQQ = y + 13*x + 16
# assert lPQ(P.xy()[0],P.xy()[1]) == lPQ(Q.xy()[0],Q.xy()[1]) == 0
# assert lPP(P.xy()[0],P.xy()[1]) == lPP((-2*P).xy()[0],(-2*P).xy()[1]) == 0
# assert lQQ(Q.xy()[0],Q.xy()[1]) == lQQ((-2*Q).xy()[0],(-2*Q).xy()[1]) == 0
# D1 = Divisor([2,1],[R,S])
# D2 = Divisor([3,-3],[R,S])
# D3 = Divisor([1,1,-2],[R,S,E(0)])
# assert evaluate_function_on_divisor(lPQ,D1) == 122
# assert evaluate_function_on_divisor(lPP,D2) == 53
