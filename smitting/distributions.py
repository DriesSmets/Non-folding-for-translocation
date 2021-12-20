from symfit.distributions import Gaussian as _gaussian
from symfit import Parameter
import numpy as np
from itertools import count


class BaseDistribution(object):
    pass


class Lorentzian(BaseDistribution):
    _ids = count(0)

    def __init__(self):
        i = next(self._ids)
        self.mu_x = Parameter('mu_x_{}'.format(i))
        self.gamma = Parameter('gamma_{}'.format(i))
        _parameters = [self.mu_x, self.gamma]
        self.parameters = {p.name: p for p in _parameters}

    def make_model(self, x):
        return _lorentzian(x, self.mu_x, self.gamma)


class Gaussian(BaseDistribution):
    _ids = count(0)

    def __init__(self):
        i = next(self._ids)
        self.mu_x = Parameter('mu_x_{}'.format(i))
        self.sigma = Parameter('sigma_{}'.format(i))
        _parameters = [self.mu_x, self.sigma]
        self.parameters = {p.name: p for p in _parameters}

    def make_model(self, x):
        return _gaussian(x, self.mu_x, self.sigma)


def _lorentzian(x, mu_x, gamma):
    return (1/(2*np.pi))*(gamma / ( (x-mu_x)**2 +(0.5*gamma)**2))