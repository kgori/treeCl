__author__ = 'kgori'

import pllpy

class PLLException(Exception):
    pass


def create_instance(alignment, partitions, tree, threads=1, rns=int("0xCA55E77E", 16)):
    try:
        with open(alignment):
            pass
    except IOError as exc:
        raise exc

    if tree in ['random', 'parsimony']:
        if tree == 'random':
            instance = pllpy.pll(alignment, partitions, False, threads, rns)
        else:
            instance = pllpy.pll(alignment, partitions, True, threads, rns)
    else:
        instance = pllpy.pll(alignment, partitions, tree, threads, rns)
    return instance


def set_partition_model_parameters(instance, partition, alpha, freqs, rates, empirical_freqs, equal_freqs):
    """
    Sets parameter values for a specific partition.

    :param instance: PLL instance being modified
    :param partition: Number of the partition having its parameters set
    :param alpha: Alpha parameter of the 4-category discrete gamma rates distribution
    :param freqs: Equilibrium frequencies of states (4 for DNA, 20 for aa)
    :param rates: Relative substitution rate parameters - values of the upper triangle of 4x4 matrix,
                  so 6 numbers in all. The sixth value must be 1.0. Assume matrix is in "acgt" order.
                  Only applies to DNA data; protein models all use empirical rates.
    :param empirical_freqs: Use empirical estimates for state frequencies. Overwrites 'freqs'.
    :param equal_freqs: Set all state frequencies to 1/num_states
    :return: void
    """
    if empirical_freqs:
        freqs = instance.get_empirical_frequencies()[partition]
    elif equal_freqs:
        if instance.is_dna(partition):
            freqs = [0.25] * 4
        else:
            freqs = [0.05] * 20
    if alpha is not None:
        instance.set_alpha(alpha, partition, True)
    if freqs is not None:
        instance.set_frequencies(freqs, partition, True)
    if rates is not None:
        instance.set_rates(rates, partition, True)


def set_params_from_dict(instance, model):
    """
    Sets parameters of pll instance according to dict
    :param instance: pll instance
    :param model: dict describing pll model parameters
    :return:
    """
    p_info = model['partitions']
    for i in range(instance.get_number_of_partitions()):
        alpha = p_info[i].get('alpha')
        freqs = p_info[i].get('frequencies')
        rates = p_info[i].get('rates')
        set_partition_model_parameters(instance, i, alpha, freqs, rates, False, False)
    return instance


def pll_to_dict(instance):
    """
    Summarises parameter values from PLL instance and writes their values
    to disk in a json format file

    :param instance: PLL instance being summarised
    :param json_file: Either a filepath or a file-like stream (e.g. sys.stdout)
    :return: void
    """
    model = {'ml_tree': instance.get_tree(), 'likelihood': instance.get_likelihood(), 'partitions': {}}
    for i in range(instance.get_number_of_partitions()):
        data = {'alpha': instance.get_alpha(i), 'frequencies': instance.get_frequencies_vector(i)}
        if instance.is_dna(i):
            data['rates'] = instance.get_rates_vector(i)
        data['model'] = instance.get_model_name(i)
        data['name'] = instance.get_partition_name(i)
        model['partitions'][i] = data
    return model
