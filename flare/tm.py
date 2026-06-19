"""
duplex melting temperature with nearest neighbor thermodynamics
"""
import math
from .params import R, NN_SETS, COMPLEMENT
from . import salt as saltmod


def _revcomp(seq):
    return ''.join(COMPLEMENT[b] for b in reversed(seq))


def nn_enthalpy_entropy(seq, nn_set):
    """
    sum nearest neighbor dH and dS for a perfect duplex plus the end terms
    """
    seq = seq.upper()
    nn = nn_set['nn']
    dH = dS = 0.0
    for i in range(len(seq) - 1):
        h, s = nn[seq[i:i + 2]]
        dH += h
        dS += s
    for end in (seq[0], seq[-1]):
        if end in 'GC':
            dH += nn_set['init_gc'][0]
            dS += nn_set['init_gc'][1]
        else:
            dH += nn_set['init_at'][0]
            dS += nn_set['init_at'][1]
    if seq == _revcomp(seq):
        dH += nn_set['symmetry'][0]
        dS += nn_set['symmetry'][1]
    return dH, dS


def tm(seq, ct_molar, na_total_molar, config):
    """
    melting temperature in celsius for the oligo against its complement
    """
    seq = seq.upper()
    x = 1 if seq == _revcomp(seq) else config.get('x_factor', 4)
    nn_set = NN_SETS[config['nn_params']]
    dH, dS = nn_enthalpy_entropy(seq, nn_set)
    dS_salt = saltmod.salt_santalucia_1998_dS(dS, na_total_molar, len(seq))
    tm_k = (dH * 1000.0) / (dS_salt + R * math.log(ct_molar / x))
    return tm_k - 273.15
