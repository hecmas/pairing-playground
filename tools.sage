def log2(x):
    if x == 0:
        return 0

    r = 1
    while x > 1:
        x = x // 2
        r += 1

    return r

def embedding_degree(q,r):    
    k = 1
    while (q**k - 1) % r != 0:
        k += 1

    return k

# Find line y = mx + c passing through two points P and Q
# or vertical line y = x0 if Q = -P
# and evaluate it at a point T
def line(P, Q, T, E):
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