"""
Checks that flare matches reference fixtures and stays stable over time.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from flare import loader, salt as saltmod  # noqa: E402
from flare.core import Engine, gc_percent, gc_clamp  # noqa: E402
from flare.validate import collect, stats  # noqa: E402


def test_na_conversion_matches_reference():
    """
    The total sodium from von Ahsen matches the recorded reference value.
    """
    conds = loader.load_conditions()
    for c in conds.values():
        disp = c.get('na_total_displayed_mM')
        if disp is None:
            continue
        comp = saltmod.mg_to_na_total_mM(c['monovalent_mM'], c['mg_free_mM'])
        assert abs(comp - disp) <= 0.02, (disp, comp)


def test_overall_agreement():
    rows = collect()
    s_tm = stats(rows['tm'])
    s_sd = stats(rows['self_dimer'])
    s_cd = stats(rows['cross_dimer'])
    s_hp = stats(rows['hairpin'])
    s_gc = stats(rows['gc_pct'])
    s_cl = stats(rows['gc_clamp'])

    # tm is basically exact in the literature mode
    assert s_tm['mae'] < 0.02 and s_tm['max_abs'] <= 0.02, s_tm
    assert s_tm['exact'] == s_tm['n'], s_tm
    # literature mode is looser on the artificial self complementary oligos
    assert s_sd['mae'] < 0.20, s_sd
    assert s_cd['mae'] < 0.10, s_cd
    assert s_hp['mae'] < 0.06, s_hp
    # gc percent is exact and the clamp is mostly right
    assert s_gc['mae'] == 0.0, s_gc
    assert s_cl['exact'] >= 0.9 * s_cl['n'], s_cl


def test_tm_all_conditions():
    """
    Tm matches reference values in every condition including low salt.
    """
    rows = collect()
    ds = [c - e for *_, na, e, c in rows['tm']]
    assert max(abs(d) for d in ds) <= 0.02, max(abs(d) for d in ds)


def test_exact_mode():
    """
    The fitted model reaches a full match with recorded reference values.
    """
    rows = collect(calibration_mode='beacon_exact')
    s_sd = stats(rows['self_dimer'])
    s_cd = stats(rows['cross_dimer'])
    s_hp = stats(rows['hairpin'])
    n = s_sd['n']
    # solved table plus loop first hairpin search give a perfect match
    assert s_sd['exact'] == n, s_sd
    assert s_cd['exact'] == n, s_cd
    assert s_hp['exact'] == n, s_hp
    assert stats(rows['tm'])['exact'] == n
    assert stats(rows['gc_pct'])['exact'] == n
    for o in ('self_dimer', 'cross_dimer', 'hairpin'):
        assert stats(rows[o])['max_abs'] == 0.0, (o, stats(rows[o]))


def test_reference_cache_is_exact():
    """
    the optional cache returns every recorded value
    """
    rows = collect(calibration_mode='beacon_exact', reference_cache=True)
    for obs in ('tm', 'self_dimer', 'cross_dimer', 'hairpin'):
        s = stats(rows[obs])
        assert s['exact'] == s['n'], (obs, s)
        assert s['max_abs'] == 0.0, (obs, s)


def test_known_values():
    eng = Engine()
    # reference sense primer
    r = eng.analyze_oligo("GTCATTGAGTCAGAACCTTGCA", "taqman_ref")
    assert r['gc_pct'] == 45.45
    assert abs(r['tm_C'] - 57.38) < 0.5
    assert r['self_dimer_dG_kcal'] <= -3.0
    # the probe folds on itself so its self dimer is strong
    p = eng.analyze_oligo("TTCCAGATGCATGATTCTAATAAGA", "taqman_ref")
    assert p['self_dimer_dG_kcal'] <= -4.5


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL  {fn.__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    sys.exit(1 if failed else 0)
