"""
build the golden snapshot of every flare value, run once then leave it frozen
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from flare import loader            # noqa: E402
from flare.core import Engine       # noqa: E402

GOLDEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'golden_engine.json')

# arbitrary oligos that are not built for any structure, they lock general behaviour
ARBITRARY = [
    "ACGTACGTACGTACGT", "GGGGGAAAAATTTTTCCCCC", "ATCGGCTAGCTAGCTAGCTA",
    "TTTTTTTTTTTTTTTTTTTT", "GCATGCATGCATGCAT", "AACCGGTTAACCGGTT",
    "CAGTCAGTCAGTCAGTCAGT", "TGACATGCATGCTAGCTAGC", "CCCCCCGGGGGG",
    "ATATATGCGCGCATATAT", "GTACGTACGTACGTACGTAC", "AGAGAGAGAGAGAGAG",
    "CTCGAGCTCGAGCTCGAG", "TACGGCTAACGTTAGCCGTA", "AAAGGGCCCTTTAAAGGG",
]
ARB_CONDITIONS = "taqman_ref"


def engines():
    return {
        'literature': Engine(calibration_mode='literature'),
        'calibrated': Engine(calibration_mode='calibrated'),
    }


def build():
    """
    run both modes over every oligo and the arbitrary battery
    """
    eng = engines()
    out = {'na_total_mM': {}, 'oligos': [], 'cross_dimers': [], 'arbitrary': [],
           'arbitrary_cross': []}

    conds = loader.load_conditions()
    for cid in sorted(conds):
        out['na_total_mM'][cid] = round(eng['literature'].na_total_mM(cid), 4)

    oligos = loader.load_oligos()
    for assay in loader.load_assays():
        cond = assay['conditions']
        for o in assay['oligos']:
            seq = oligos[o['id']]
            rec = {'id': o['id'], 'seq': seq, 'conditions': cond}
            for mode, e in eng.items():
                r = e.analyze_oligo(seq, cond)
                rec[mode] = {k: r[k] for k in ('tm_C', 'gc_pct', 'gc_clamp',
                                               'self_dimer_dG_kcal', 'hairpin_dG_kcal')}
            out['oligos'].append(rec)
        for cd in assay.get('cross_dimers', []):
            a, b = cd['pair']
            rec = {'assay': assay['id'], 'pair': [a, b]}
            for mode, e in eng.items():
                rec[mode] = e.cross_dimer_dG(oligos[a], oligos[b])[0]
            out['cross_dimers'].append(rec)

    for seq in ARBITRARY:
        rec = {'seq': seq}
        for mode, e in eng.items():
            rec[mode] = e.analyze_oligo(seq, ARB_CONDITIONS)
        out['arbitrary'].append(rec)
    for i in range(len(ARBITRARY)):
        for j in range(i + 1, len(ARBITRARY)):
            rec = {'a': ARBITRARY[i], 'b': ARBITRARY[j]}
            for mode, e in eng.items():
                rec[mode] = e.cross_dimer_dG(ARBITRARY[i], ARBITRARY[j])[0]
            out['arbitrary_cross'].append(rec)
    return out


if __name__ == '__main__':
    with open(GOLDEN, 'w') as f:
        json.dump(build(), f, indent=1, sort_keys=True)
    print("wrote", GOLDEN)
