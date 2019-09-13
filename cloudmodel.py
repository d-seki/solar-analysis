from scipy.optimize import minimize, basinhopping
import numpy as np


class CloudModel:
    def __init__(self, Cl, bg2, intl_prms):
        self.bg2 = bg2
        if bg2.shape[0] == 73:
            self.lobs = np.linspace(-9.,9.,bg2.shape[0])
        elif bg2.shape[0] == 25:
            self.lobs = np.linspace(-3.,3.,bg2.shape[0])
        self.Cl = Cl
        self.intl_prms = intl_prms
        self.co = 2
        #self.intl_prms = np.array([0.2, 0.2, 0.3, 0.])
        # l0/c*np.sqrt(2*kB*T/MH) = 0.28

    def cntrst_func(self, params):
        S = params[0]
        tau = params[1]
        lwid = params[2]
        ldop = params[3]
        taul = tau*np.exp(-(((self.lobs-ldop)/lwid)**2))
        return (S/self.bg2 - 1.)*(1 - np.exp(-taul))

    def residual(self, params):
        res_array = self.Cl - self.cntrst_func(params)
        return (res_array*res_array).sum()

    def compute(self):
        kwargs = {'method': 'Powell'}
        sol = minimize(self.residual, self.intl_prms, method='Powell')
        #sol = basinhopping(self.residual, self.intl_prms,minimizer_kwargs=kwargs)
        return sol

#CPU times: user 1min 40s, sys: 392 ms, total: 1min 40s
#Wall time: 1min 42s
