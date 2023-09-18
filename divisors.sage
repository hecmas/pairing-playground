
class Divisor():
    def __init__(self, multiplicites, points):
        assert len(multiplicites) == len(points)
        for m in multiplicites:
            assert m in ZZ
        # This is not necessary, assuming that every point P in points is defined as E(P)
        # for P in points:
        #     assert P in E

        if len(points) >= 1:
            self.E = points[0].curve()

            for i in range(1,len(points)):
                if points[i].curve() != self.E:
                    raise ValueError("Points must be defined over the same curve")
        else:
            self.E = None

        for i in range(len(multiplicites)-1,-1,-1):
            if multiplicites[i] == 0:
                multiplicites.pop(i)
                points.pop(i)

        self.multiplicites = multiplicites
        self.points = points
        self.degree = sum(multiplicites)
        self.support = points

        self.repr = {}
        for i in range(len(multiplicites)):
            self.repr[points[i]] = multiplicites[i]

    def __repr__(self):
        if len(self.multiplicites) == 0:
            return "0"

        result = ""
        for i in range(len(self.multiplicites)):
            if self.multiplicites[i] > 0:
                if i > 0:
                    result += " + "
                result += str(self.multiplicites[i])
            else:
                result += " - "  + str(self.multiplicites[i])[1:]

            result += "·" + str(self.points[i].xy())

        return result

    def __add__(self, other):
        multiplicites = []
        points = []
        for key in self.repr.keys():
            if key in other.repr.keys():
                keysum = self.repr[key] + other.repr[key]
                if keysum != 0:
                    multiplicites.append(keysum)
                    points.append(key)
            else:
                multiplicites.append(self.repr[key])
                points.append(key)

        for key in other.repr.keys():
            if key not in self.repr.keys():
                multiplicites.append(other.repr[key])
                points.append(key)

        return Divisor(multiplicites, points)

    def __neg__(self):
        negative = [-m for m in self.multiplicites]
        return Divisor(negative, self.points)

    def __sub__(self, other):
        return self + (-other)

    def __eq__(self, other):
        if len(self.multiplicites) != len(other.multiplicites):
            return False

        for m in self.multiplicites:
            if m not in other.multiplicites:
                return False
            
        for P in self.points:
            if P not in other.points:
                return False

        return True

    def __ne__(self, other):
        return not self == other

    def is_principal(self):
        if self.is_zero():
            return True

        sum = self.E(0)
        for i in range(len(self.multiplicites)):
            sum += self.multiplicites[i] * self.points[i]

        if sum == self.E(0) and self.degree == 0:
            return True
        else:
            return False
        
    def equivalent(self, other):
        D = self - other
        return D.is_principal()

    def is_zero(self):
        return self == Divisor.zero()

    @classmethod
    def zero(cls):
        return Divisor([], [])


def evaluate_function_on_divisor(f, D):
    # TODO: Add a check for disjoint supports
    #       Implement f --> (f), it seems Sage does not support this

    result = 1
    for i in range(len(D.multiplicites)):
        if D.points[i].is_zero():
            # homogenize the function and evaluate at the origin in projective space
            F = f.numerator().homogenize() / f.denominator().homogenize()
            result *= F(0, 1, 0) ^ D.multiplicites[i]
        else:
            result *= f(D.points[i].xy()) ^ D.multiplicites[i]

    return result
