import pyprelude.FPToolBox as fp
from math import factorial
from sympy import *
from sympy.abc import a, b, x

# symbolic gate function
#eta = lambda a, b, x: Piecewise((1, (x >= a) & (x <= b)), (0, True))

class eta(Function):
    nargs = 3

    @classmethod
    def eval(cls, a, b, x):
        if a.is_Number and b.is_Number and x.is_Number:
            zero = numbers.Zero()
            if x >= a and x <=b:
                return 1 + zero
            else:
                return 0 + zero

# symbolic theshold function
#theta = lambda a, b, x: Piecewise((a, x < a), (b, x > b), (x, True))

class theta(Function):
    nargs = 3

    @classmethod
    def eval(cls, a, b, x):
        if a.is_Number and b.is_Number and x.is_Number:
            if x < a:
                return a
            elif x >b:
                return b
            else:
                return x

# get unique representation of an edge
def e_repr(e):
    if e[0]<e[1]:
        return e
    else:
        return (e[1], e[0])

# return ordered set for B_curl(a,b)
# where a is the dimension of the tuple, and k is the sum
def B_curl(n, k):
    tmp = B_curl_h(n, k)
    return list(map(tuple, tmp))

def B_curl_h(n, k):
    if n == 1:
        return [[k]]

    if k == 0:
        return fp.lmap(lambda xs: addhd(0, xs), B_curl_h(n - 1, 0))
    
    # for each possible value, get one
    res = []
    for i in range(k + 1):
        res += fp.lmap(lambda xs: addhd(i, xs), B_curl_h(n - 1, k - i))

    return res

# the number of elements in B_curl
def B_curl_num(n, k):
    if n == 1:
        return 1

    if k == 0:
        return 1
    
    res = 0
    for i in range(k + 1):
        res += B_curl_num(n - 1, k - i)

    return res

# return ordered set for A_curl(a,b)
def A_curl(n, k):
    res = []
    for k0 in range(k + 1):
        res += B_curl(n, k0)

    return res

# return ordered set for A_curl(a,b)
def A_curl_num(n, k):
    res = 0
    for k0 in range(k + 1):
        res += B_curl_num(n, k0)

    return res

# function to fast add head to a list
def addhd(hd, xs):
    xs.insert(0, hd)
    return xs

#def ncr(n, r):
#    r = min(r, n-r)
#    numer = reduce(op.mul, range(n, n-r, -1), 1)
#    denom = reduce(op.mul, range(1, r+1), 1)
#    return numer / denom

# calc n choose multiindex rs

def ncrs(n, rs):
    res = factorial(n)
    r_sum = 0
    for r in rs:
        res = res / factorial(r)
        r_sum = r_sum + r

    return res / factorial(n - r_sum)

# search symetric tuples
#def get(mydict, key):
#    if key[0] > key[1]:
#        return mydict[key[1], key[0], key[2], key[3]]
#    else:
#        return mydict[key]

# convert to sympy expression
def toSym(val):
    zero = numbers.Zero()
    return val + zero

# rationalize input
def rat(val, rational):
    if rational is False:
        return float(val)
    else:
        return Rational(val)


if __name__ == '__main__':
    from sympy.plotting import plot
    
    plot(eta(1, oo, x), (x, -1, 5))
