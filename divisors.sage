
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
        if len(self.multiplicites) == 0:
            return "0"

        result = ""
        for i in range(len(self.multiplicites)-1):
            result += str(self.multiplicites[i]) + " " + str(self.points[i].xy()) + " + "

        result += str(self.multiplicites[-1]) + " " + str(self.points[-1].xy())
        
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
        negative = [-m for m in self.multiplicites]
        return Divisor(negative, self.points)

    def __sub__(self, other):
        return self + (-other)
    
    def __eq__(self, other):
        return self.sum == other.sum
    
    def __ne__(self, other):
        return not self == other

    def is_principal(self, E):
        addition = E(0)
        for i in range(len(self.multiplicites)):
            addition += self.multiplicites[i] * self.points[i]

        if addition == E(0) and self.degree == 0:
            return True
        else:
            return False

    def is_zero(self):
        return self == Divisor.zero()
    
    @classmethod
    def zero(cls):
        return Divisor([], [])

def evaluate_function_on_divisor(f,D):
    # TODO: Add a check for disjoint supports
   
    result = 1
    for i in range(len(D.multiplicites)):
        if D.points[i].is_zero():
            # homogenize the function and evaluate at the origin in projective space
            F = f.numerator().homogenize() / f.denominator().homogenize()
            result *= F(0,1,0) ^ D.multiplicites[i]
        else:
            result *= f(D.points[i].xy()) ^ D.multiplicites[i]

    return result
