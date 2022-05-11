import numpy as np
from scipy.optimize import newton


def viscousSublayer(yPlus):
    return yPlus


def logLowRegion(yPlus, kappa=0.41, B=5.2):
    yPlus = np.where(np.abs(yPlus)>1e-30, yPlus, 1e-30)
    return np.log(yPlus) / kappa + B


def wakeContribution(y, delta=1., kappa=0.41, Pi=0.18):
    return 2 * np.sin(np.pi/2 * y/delta)**2 / kappa * Pi


def wallLow(yPlus, kappa=0.41, B=5.2, wake=False, delta=1., Pi=0.18):
    func = lambda yPlus: viscousSublayer(yPlus) - logLowRegion(yPlus, kappa, B)
    yPlus_ = newton(func=func, x0=11)
    if wake==True:
        return np.where(yPlus<yPlus_, viscousSublayer(yPlus), logLowRegion(yPlus, kappa, B)) + wakeContribution(y, delta, kappa, Pi)
    else:
        return np.where(yPlus<yPlus_, viscousSublayer(yPlus), logLowRegion(yPlus, kappa, B))
    

if __name__=="__main__":
    import matplotlib.pyplot as plt

    Ret = 400.
    y = np.linspace(0, 1, int(Ret*2))
    yPlus = y * Ret

    plt.figure()
    plt.plot(yPlus, wallLow(yPlus))
    plt.plot(yPlus, wallLow(yPlus, wake=True))
    plt.grid()
    plt.xlim(1, Ret)
    plt.ylim(0, 30)
    plt.xscale('log')
    plt.show()
