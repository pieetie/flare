"""
thermodynamic tables for tm and dimer energies
"""

R = 1.98720425864083

# santalucia 1998 nearest neighbor dH (kcal/mol) and dS (cal/mol/K)
NN_SANTALUCIA_1998 = {
    'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2),
    'AT': (-7.2, -20.4),
    'TA': (-7.2, -21.3),
    'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7),
    'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
    'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0),
    'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
    'CG': (-10.6, -27.2),
    'GC': (-9.8, -24.4),
    'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9),
}
INIT_GC_SANTALUCIA_1998 = (0.1, -2.8)
INIT_AT_SANTALUCIA_1998 = (2.3, 4.1)
SYMMETRY_SANTALUCIA_1998 = (0.0, -1.4)

NN_SETS = {
    'santalucia_2004': {
        'nn': NN_SANTALUCIA_1998,
        'init_gc': INIT_GC_SANTALUCIA_1998,
        'init_at': INIT_AT_SANTALUCIA_1998,
        'symmetry': SYMMETRY_SANTALUCIA_1998,
    },
}

# dimer dG37 calibrated fit, one value per nearest neighbor step
DG37_CALIBRATED_FIT_RAW = {
    'AA': -1.1370, 'AT': -0.9730, 'TA': -0.7130, 'CA': -1.5870, 'GT': -1.5730,
    'CT': -1.3955, 'GA': -1.4430, 'CG': -2.3500, 'GC': -2.3830, 'GG': -1.9250,
}
DG37_CALIBRATED_FIT_INIT = 1.87
DG37_CALIBRATED_FIT_AT_PENALTY = 0.1455
DG37_CALIBRATED_FIT_THRESHOLD = 0.35
DG37_CALIBRATED_FIT_HAIRPIN_NUCLEATION = 1.88
DG37_CALIBRATED_FIT_HAIRPIN_AT_PENALTY = 0.14
DG37_CALIBRATED_FIT_HAIRPIN_THRESHOLD = 0.0

COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}


def _expand_dg37(base):
    """
    fill the reverse complement steps so any dinucleotide can be looked up
    """
    out = dict(base)
    for k, v in base.items():
        rc = ''.join(COMPLEMENT[b] for b in reversed(k))
        out.setdefault(rc, v)
    return out


DG37_CALIBRATED_FIT = _expand_dg37(DG37_CALIBRATED_FIT_RAW)
