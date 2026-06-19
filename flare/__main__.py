"""
command line entry, pass the sense primer, the antisense primer and the probe
"""
import sys
from .core import Engine

CONDITIONS = "taqman_ref"


def run(sense, antisense, probe):
    """
    print the values for the three oligos and their cross dimers
    """
    eng = Engine(calibration_mode="calibrated")
    seqs = {"sense": sense, "antisense": antisense, "probe": probe}
    print(f"{'oligo':10}{'len':>4}{'Tm':>8}{'GC%':>7}{'clamp':>6}{'self':>7}{'hairpin':>9}")
    for name, s in seqs.items():
        r = eng.analyze_oligo(s, CONDITIONS)
        print(f"{name:10}{r['length']:>4}{r['tm_C']:>8.2f}{r['gc_pct']:>7.2f}"
              f"{r['gc_clamp']:>6}{r['self_dimer_dG_kcal']:>7.1f}{r['hairpin_dG_kcal']:>9.1f}")
    print()
    for x, y in (("sense", "antisense"), ("sense", "probe"), ("antisense", "probe")):
        dg, _ = eng.cross_dimer_dG(seqs[x], seqs[y])
        print(f"cross {x} {y} {dg:+.1f}")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 3:
        print("usage  python -m flare SENSE ANTISENSE PROBE")
        return 1
    run(*argv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
