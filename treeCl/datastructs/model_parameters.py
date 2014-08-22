from abc import ABCmeta
import json

class ModelParameters:
    __metaclass__ = ABCmeta
    def __init__(self, name, model):
        self.name = name
        self.model = model

    def __repr__(self):
        return super(ModelParameters, self).__repr__()

    def __str__(self):
        return super(ModelParameters, self).__str__()

    def set_name(self, name):
        self._name = name

    def get_name(self):
        return self._name

    name = abstractproperty(get_name, set_name)

    def set_model(self, model):
        self._model = model

    def get_model(self):
        return self._model

    model = abstractproperty(get_model, set_model)

    def set_frequencies(self, freqs):
        self._freqs = freqs

    def get_frequencies(self):
        return self._freqs

    frequencies = abstractproperty(get_frequencies, set_frequencies)

    def set_rates(self, rates):
        self._rates = rates

    def get_rates(self):
        return self._rates

    rates = abstractproperty(get_rates, set_rates)

    def set_alpha(self, alpha):
        self._alpha = alpha

    def get_alpha(self):
        return self._alpha

    alpha = abstractproperty(get_alpha, set_alpha)

    @abstractmethod
    def get_json_string(self):
        pass