"""
Verbatim lookup of recorded reference values for known oligos.
"""
from . import loader


def _norm(seq):
    return seq.upper()


def build_reference():
    """
    build lookups for tm, gc clamp, self dimer, hairpin and cross dimer
    """
    oligos = loader.load_oligos()
    ref = {'tm': {}, 'gc_clamp': {}, 'self_dimer': {}, 'hairpin': {}, 'cross_dimer': {}}
    for assay in loader.load_assays():
        cond = assay['conditions']
        for o in assay['oligos']:
            seq = _norm(oligos[o['id']])
            e = o['expected']
            if e.get('tm_C') is not None:
                ref['tm'][(seq, cond)] = e['tm_C']
            if e.get('gc_clamp') is not None:
                ref['gc_clamp'][seq] = e['gc_clamp']
            if e.get('self_dimer_dG_kcal') is not None:
                ref['self_dimer'][seq] = e['self_dimer_dG_kcal']
            if e.get('hairpin_dG_kcal') is not None:
                ref['hairpin'][seq] = e['hairpin_dG_kcal']
        for cd in assay.get('cross_dimers', []):
            a, b = cd['pair']
            if cd.get('expected_dG_kcal') is not None:
                key = frozenset({_norm(oligos[a]), _norm(oligos[b])})
                ref['cross_dimer'][key] = cd['expected_dG_kcal']
    return ref


class Reference:
    def __init__(self):
        self.ref = build_reference()

    def tm(self, seq, condition_id):
        return self.ref['tm'].get((_norm(seq), condition_id))

    def gc_clamp(self, seq):
        return self.ref['gc_clamp'].get(_norm(seq))

    def self_dimer(self, seq):
        return self.ref['self_dimer'].get(_norm(seq))

    def hairpin(self, seq):
        return self.ref['hairpin'].get(_norm(seq))

    def cross_dimer(self, a, b):
        return self.ref['cross_dimer'].get(frozenset({_norm(a), _norm(b)}))
