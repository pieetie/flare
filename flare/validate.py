"""
Compare engine output against recorded reference values.
"""
import math
from . import loader
from .core import Engine, gc_percent, gc_clamp

OBSERVABLES = ['tm', 'gc_pct', 'gc_clamp', 'self_dimer', 'cross_dimer', 'hairpin']


def collect(calibration_mode='literature'):
    """
    gather expected and computed pairs for every observable
    """
    eng = Engine(calibration_mode=calibration_mode)
    oligos = loader.load_oligos()
    rows = {k: [] for k in OBSERVABLES}
    for assay in loader.load_assays():
        cond = assay['conditions']
        na = round(eng.na_total_mM(cond), 2)
        for o in assay['oligos']:
            seq = oligos[o['id']]
            e = o['expected']
            tag = (assay['id'], o['id'], na)
            if e.get('tm_C') is not None:
                rows['tm'].append((*tag, e['tm_C'], eng.tm(seq, cond)))
            if e.get('gc_pct') is not None:
                rows['gc_pct'].append((*tag, e['gc_pct'], gc_percent(seq)))
            if e.get('gc_clamp') is not None:
                rows['gc_clamp'].append((*tag, e['gc_clamp'], gc_clamp(seq)))
            if e.get('self_dimer_dG_kcal') is not None:
                rows['self_dimer'].append((*tag, e['self_dimer_dG_kcal'], eng.self_dimer_dG(seq)[0]))
            if e.get('hairpin_dG_kcal') is not None:
                rows['hairpin'].append((*tag, e['hairpin_dG_kcal'], eng.hairpin_dG(seq)[0]))
        for cd in assay.get('cross_dimers', []):
            a, b = cd['pair']
            if cd.get('expected_dG_kcal') is not None:
                rows['cross_dimer'].append((assay['id'], f"{a}x{b}", na,
                                            cd['expected_dG_kcal'], eng.cross_dimer_dG(oligos[a], oligos[b])[0]))
    return rows


def stats(rows):
    """
    error summary for one observable
    """
    ds = [comp - exp for *_, exp, comp in rows]
    n = len(ds)
    if n == 0:
        return None
    return dict(n=n, mae=sum(abs(d) for d in ds) / n,
                rmse=math.sqrt(sum(d * d for d in ds) / n),
                bias=sum(ds) / n, max_abs=max(abs(d) for d in ds),
                exact=sum(1 for d in ds if abs(d) < 0.05))


def main():
    for mode in ('literature', 'calibrated'):
        rows = collect(mode)
        print(mode)
        for obs in OBSERVABLES:
            s = stats(rows[obs])
            print(f"  {obs:12} exact {s['exact']:3}/{s['n']:<3} mae {s['mae']:.3f} max {s['max_abs']:.2f}")


if __name__ == '__main__':
    main()
