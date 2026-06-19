"""
locks the public api, the tuned constants, determinism and the edge cases
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from flare import params, salt as saltmod  # noqa: E402
from flare.core import Engine, gc_percent, gc_clamp  # noqa: E402


# the calibrated constants must not move
def test_constants_na_and_tm():
    assert saltmod.VON_AHSEN_K == 126.49
    cfg = Engine().cfg
    assert cfg['nn_params'] == 'santalucia_2004'
    assert cfg['salt_correction'] == 'santalucia_1998'
    assert cfg['x_factor'] == 1


def test_constants_fit_dimer_table():
    expected = {
        'AA': -1.1370, 'AT': -0.9730, 'TA': -0.7130, 'CA': -1.5870, 'GT': -1.5730,
        'CT': -1.3955, 'GA': -1.4430, 'CG': -2.3500, 'GC': -2.3830, 'GG': -1.9250,
    }
    for k, v in expected.items():
        assert params.DG37_CALIBRATED_FIT_RAW[k] == v, (k, params.DG37_CALIBRATED_FIT_RAW[k])
    assert params.DG37_CALIBRATED_FIT_INIT == 1.87
    assert params.DG37_CALIBRATED_FIT_AT_PENALTY == 0.1455
    assert params.DG37_CALIBRATED_FIT_THRESHOLD == 0.35
    assert params.DG37_CALIBRATED_FIT_HAIRPIN_NUCLEATION == 1.88
    assert params.DG37_CALIBRATED_FIT_HAIRPIN_AT_PENALTY == 0.14


def test_santalucia_1998_nn_table():
    # the published nn values that drive tm
    nn = params.NN_SANTALUCIA_1998
    assert nn['AA'] == (-7.9, -22.2)
    assert nn['CG'] == (-10.6, -27.2)
    assert nn['GC'] == (-9.8, -24.4)


# the public methods and helpers must stay callable
def test_engine_api_surface():
    for mode in ('literature', 'calibrated'):
        e = Engine(calibration_mode=mode)
        for m in ('tm', 'na_total_mM', 'self_dimer_dG', 'cross_dimer_dG',
                  'hairpin_dG', 'analyze_oligo'):
            assert callable(getattr(e, m)), (mode, m)
    r = Engine().analyze_oligo("GTCATTGAGTCAGAACCTTGCA", "taqman_ref")
    assert set(r) >= {'sequence', 'length', 'tm_C', 'gc_pct', 'gc_clamp',
                      'self_dimer_dG_kcal', 'hairpin_dG_kcal'}
    assert gc_percent("GGCC") == 100.0
    assert isinstance(gc_clamp("GTCATTGAGTCAGAACCTTGCA"), int)


def test_modes_exist():
    Engine(calibration_mode='literature')
    Engine(calibration_mode='calibrated')
    try:
        Engine(calibration_mode='nonsense')
        assert False, "should reject unknown mode"
    except ValueError:
        pass


# same input gives the same output and short oligos do not crash
def test_determinism():
    e = Engine(calibration_mode='calibrated')
    a = e.analyze_oligo("ACGTACGTACGTACGT", "taqman_ref")
    b = e.analyze_oligo("ACGTACGTACGTACGT", "taqman_ref")
    assert a == b


def test_case_insensitive():
    e = Engine()
    assert e.tm("gtcattgagtcagaaccttgca", "taqman_ref") == e.tm("GTCATTGAGTCAGAACCTTGCA", "taqman_ref")


def test_edge_cases_no_crash():
    e = Engine(calibration_mode='calibrated')
    for s in ("AT", "GC", "AAAAA", "GGGGGGGG", "ATGCATGCAT", "A"):
        r = e.analyze_oligo(s, "taqman_ref")
        for key in ('self_dimer_dG_kcal', 'hairpin_dG_kcal'):
            assert isinstance(r[key], float)


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    bad = 0
    for fn in fns:
        try:
            fn(); print("PASS ", fn.__name__)
        except Exception as e:
            bad += 1; print("FAIL ", fn.__name__, ":", e)
    print(f"\n{len(fns)-bad}/{len(fns)} passed")
    sys.exit(1 if bad else 0)
