import numpy as np

def _preprocess_inputs(x, weights):
    """
    Coerce inputs into compatible format
    """
    if weights is None:
        w_arr = np.ones(len(x))
    else:
        w_arr = np.array(weights)
    x_arr = np.array(x)
    if x_arr.ndim == 2:
        if w_arr.ndim == 1:
            w_arr = w_arr[:, np.newaxis]
    return x_arr, w_arr

def amean(x, weights=None):
    """
    Return the weighted arithmetic mean of x
    """
    w_arr, x_arr = _preprocess_inputs(x, weights)
    return (w_arr*x_arr).sum(axis=0) / w_arr.sum(axis=0)

def gmean(x, weights=None):
    """
    Return the weighted geometric mean of x
    """
    w_arr, x_arr = _preprocess_inputs(x, weights)
    return np.exp((w_arr*np.log(x_arr)).sum(axis=0) / w_arr.sum(axis=0))

def hmean(x, weights=None):
    """
    Return the weighted harmonic mean of x
    """
    w_arr, x_arr = _preprocess_inputs(x, weights)
    return w_arr.sum(axis=0) / (w_arr/x_arr).sum(axis=0)
