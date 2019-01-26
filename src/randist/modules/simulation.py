from sympy.stats import FiniteRV, Uniform, sample

def min_dist(ev1, ev2, g):
    e = strfy(ev1[0])
    p = ev1[1]
    f = strfy(ev2[0])
    q = ev2[1]
    # on the same edge
    le = g.l[e]
    lf = g.l[f]
    if e == f:
        duv = g.d[e[0]][e[1]]
        d0 = abs(p - q) * le
        d1 = (1 - abs(p - q)) * le + duv
        res = min(d0, d1)
    else:
        d00 = g.d[e[0]][f[0]] + le * p + lf * q
        d01 = g.d[e[0]][f[1]] + le * p + lf * (1 - q)
        d10 = g.d[e[1]][f[0]] + le * (1 - p) + lf * q
        d11 = g.d[e[1]][f[1]] + le * (1 - p) + lf * (1 - q)
        res = min(min(min(d00, d01), d10), d11)

    return res

def strfy(tl):
    return (str(tl[0]), str(tl[1]))

class events:
    def __init__(self, g):
        self.g = g
        x_den = {}
        y_den = {}
        for e in g.edges():
            x_den[e] = g.gx[e]
            y_den[e] = g.gy[e]


        self.X = FiniteRV('X', x_den)
        self.Y = FiniteRV('Y', y_den)
        self.P = Uniform('P', 0, 1)
        self.Q = Uniform('Q', 0, 1)

    def x_sample(self):
        x = sample(self.X)
        p = sample(self.P)
        return(x, p)

    def y_sample(self):
        y = sample(self.Y)
        q = sample(self.Q)
        return(y, q)

    def d_sample(self):
        ev1 = self.x_sample()
        ev2 = self.y_sample()
        return min_dist(ev1, ev2, self.g)

    def E_d(self, f=lambda x: x, nsample= 10000):
        res = 0
        for i in range(nsample):
            if i % 1000 == 0:
                print(i)
            res = i / (i + 1.0) * res + f(self.d_sample()) / (i + 1.0)

        return res
            
    




        
