import numpy as np
import matplotlib.pyplot as plt
from .pyprelude import FPToolBox as fp

def plot1d(func, lb, ub, step):
    x = np.arange(lb, ub, step)
    y = fp.lmap(func, x)
    plt.figure(1)
    plt.minorticks_off()
    plt.plot(x, y, 'k', color='b')
    plt.show()
