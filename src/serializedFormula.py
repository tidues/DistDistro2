from enums import Stats
import numericFuncs as nf
from sympy import *
from sympy.abc import p, q, x
from plot import plot1d

# serialized formula
class SFormula:
    def __init__(self):
        self.stat = None
        self.fs = {}
        self.N_fs = {}
        self.mods = None
        self.idx_num = 0  # dimension of indexes for self.fs
        self.val_keys = None  # the symbols for values
        self.plot_info = (None, x, 0, 0)  # plotting parameters
        self.unikey = None

    def eval(self, *params):
        # split params into keys and vals
        keys, subval, vals = self.__split_input(*params)
        keys = self.mkkey(keys)
        if keys in self.N_fs:
            myf = self.N_fs[keys]
            return self.__apply(myf, vals)
        else:
            print('This function does not exist!')
            return False

    # output the formula
    def formula(self, *params):
        # split params into keys and vals
        keys, subval, vals = self.__split_input(*params)
        keys = self.mkkey(keys)

        if keys in self.fs:
            myf = self.fs[keys]
            return self.__subs(myf, subval)
        else:
            print('This function does not exist!')
            return False

    # input f_info = (f, var, lb, ub)
    def plot(self, *params, step=0.01):
        if self.stat != Stats.MOMENT:
            if len(params) == 0:
                if self.unikey in self.fs:
                    f = self.fs[self.unikey]
                else:
                    print('This function does not exist!')
                    return False
            else:
                # split params into keys and vals
                keys, subval, vals = self.__split_input(*params)
                keys = self.mkkey(keys)
                if keys in self.fs:
                    f = self.__subs(self.fs[keys], subval)
                else:
                    print('This function does not exist!')
                    return False

            var = self.plot_info[1]
            lb = self.plot_info[2]
            ub = self.plot_info[3]

            f_lambda = lambdify(var, f, modules=self.mods)
            plot1d(f_lambda, lb, ub, step)

        else:
            print('moment function cannot be plotted.')

    def __apply(self, f, vals):
        if len(vals) == 0:
            return f
        else:
            return f(*vals)

    def mkkey(self, keys):
        if len(keys) == 1:
            return keys[0]
        else:
            return keys

    def __split_input(self, *params):
        if self.idx_num == 0:
            keys = self.unikey
            vals = params
        else:
            keys = []
            vals = []
            for i, v in enumerate(params):
                if i < self.idx_num:
                    keys.append(v)
                else:
                    vals.append(v)
            keys = tuple(keys)
        
        # make subs info
        if len(vals) == 0:
            subval = None
        elif len(vals) == 1:
            subval = [(self.val_keys[0], vals[0])]
        else:
            subval = list(zip(self.val_keys, vals))

        vals = tuple(vals)

        return (keys, subval, vals)

    # return formula with subs
    def __subs(self, formula, subval):
        if subval is None:
            return formula
        else:
            return formula.subs(subval)

