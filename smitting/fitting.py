from symfit import Parameter, Variable, parameters, Fit, MatrixSymbol, CallableModel, BlockMatrix, D, ODEModel
import numpy as np
from smitting.distributions import BaseDistribution
from smitting.interactive import Interactive
from functools import reduce
from operator import add


class BaseFit(object):
    def __init__(self):
        self.res = None
        self.model = None
        self.parameters = {}
        self.variables = {}
        self.fit = None

        self.independent_data = {}
        self.dependent_data = {}

    def make_fit(self, **kwargs):
        self.fit = Fit(self.model, **self.dependent_data, **self.independent_data, **kwargs)

    def execute(self, **kwargs):
        self.res = self.fit.execute(**kwargs)

    def set_par(self, name, value=None, min=None, max=None, fixed=None):
        par = self.parameters[name]
        for attr in ['value', 'min', 'max', 'fixed']:
            if locals()[attr] is not None:
                setattr(par, attr, locals()[attr])

    @property  #chache?
    def ans(self):
        if self.res is not None:
            return self.model(**self.independent_data, **self.res.params)
        else:
            return None


class ODEFit(BaseFit, Interactive):
    operators = ['<->', '<-', '->']

    def __init__(self, time, **kwargs):
        super(ODEFit, self).__init__()
        self.independent_data = {'t': time}
        self.dependent_data = kwargs

        self.t = Variable('t')
        self.variables['t'] = self.t

    def make_model(self, string, initial, **kwargs):
        split = string.split(' ')
        states = [s for s in split if not s in self.operators]

        model_dict = {}
        for state in states:
            v = self._get_var(state)
            key = D(v, self.t)
            components = []

            i = split.index(state)
            # look before i:
            try:
                op = split[i - 1]
                target_state = split[i - 2]  #refactor to other state
                if op in ['->', '<->']:  # Formation from target state to state
                    name = 'k_{}{}'.format(target_state, state)
                    p = self._get_param(name)
                    components.append(10**p * self._get_var(target_state))
                if op in ['<-', '<->']:  # Formation from state to target state
                    name = 'k_{}{}'.format(state, target_state)
                    p = self._get_param(name)
                    components.append(-10**p * v)
            except IndexError:
                pass

            # look after i
            try:
                op = split[i + 1]
                target_state = split[i + 2]
                if op in ['<-', '<->']:  # Formation from target state to state
                    name = 'k_{}{}'.format(target_state, state)
                    p = self._get_param(name)
                    components.append(10**p * self._get_var(target_state))
                if op in ['->', '<->']:  # Formation from state to target state
                    name = 'k_{}{}'.format(state, target_state)
                    p = self._get_param(name)
                    components.append(-10**p * v)
            except IndexError:
                pass

            model_dict[key] = reduce(add, components)

        initial_dict = {self.variables[k]: v for k, v in initial.items()}
        self.model = ODEModel(model_dict, initial_dict, **kwargs)

    def _get_param(self, name):
        try:
            return self.parameters[name]
        except KeyError:
            p = Parameter(name)
            self.parameters[name] = p
            return p

    def _get_var(self, name):
        try:
            return self.variables[name]
        except KeyError:
            v = Variable(name)
            self.variables[name] = v
            return v


class LMFit(BaseFit):
    #todo parameters name lookup object thingy
    def __init__(self, x, y):
        super(LMFit, self).__init__()

        x = x if x.ndim == 2 else np.atleast_2d(x).T

        self.independent_data = {'x': x}
        self.N_x = len(x)
        self.dependent_data = {'Y': y}
        self.N_t = len(y.T)

    def make_model(self, bases):
        """

        :param bases: dict of name: base
        :return:
        """
        A_param = Parameter('A')
        A_param.min = 0
        self.parameters = {'A': A_param}

        Y = MatrixSymbol('Y', self.N_x, self.N_t)
        x = Variable('x')

        bases_symbols = {k: MatrixSymbol(k, self.N_x, 1) for k in bases.keys()}
        A = MatrixSymbol(A_param, len(bases), self.N_t)
        B = BlockMatrix([list(bases_symbols.values())])
        model_dict = {Y: B * A}

        #todo names will clash when mixing lorentzians with gaussians
        for k, v in bases.items():
            if type(v) == type and issubclass(v, BaseDistribution):
                dist_instance = v()
                m = dist_instance.make_model(x)
                self.parameters.update(dist_instance.parameters)

                model_dict[bases_symbols[k]] = m
            else:
                _v = self._check_shape(v)
                if _v is None:
                    raise ValueError('Base {} has an invalid shape'.format(k))

                self.independent_data[k] = _v

        self.model = CallableModel(model_dict)

    def _check_shape(self, arr):
        if arr.ndim == 1 and len(arr) == self.N_x:
            return np.atleast_2d(arr).T
        elif arr.ndim == 2 and arr.shape == (self.N_x, 1):
            return arr
        else:
            return None

