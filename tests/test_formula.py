"""
The standalone single-file module must match the package exact mode.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import flare_beacon as fb                 # noqa: E402
from flare import loader                  # noqa: E402
from flare.core import Engine             # noqa: E402


def test_standalone_matches_flare():
    """
    same values from the single file and the package on every oligo
    """
    oligos = loader.load_oligos()
    conds = loader.load_conditions()
    eng = Engine(calibration_mode='beacon_exact')
    for assay in loader.load_assays():
        c = conds[assay['conditions']]
        kw = dict(dna_nM=c['nucleic_acid_nM'], monovalent_mM=c['monovalent_mM'],
                  mg_free_mM=c['mg_free_mM'])
        for o in assay['oligos']:
            s = oligos[o['id']]
            r1 = eng.analyze_oligo(s, assay['conditions'])
            r2 = fb.analyze(s, **kw)
            for k in ('tm_C', 'gc_pct', 'gc_clamp', 'self_dimer_dG_kcal', 'hairpin_dG_kcal'):
                assert r1[k] == r2[k], (o['id'], k, r1[k], r2[k])
        for cd in assay.get('cross_dimers', []):
            a, b = cd['pair']
            assert eng.cross_dimer_dG(oligos[a], oligos[b])[0] == fb.cross_dimer(oligos[a], oligos[b])


if __name__ == '__main__':
    test_standalone_matches_flare()
    print("PASS  test_standalone_matches_flare")
