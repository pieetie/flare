"""
command line entry, pass the sense primer, the antisense primer and the probe
"""
import sys
from .core import Engine
from . import render

CONDITIONS = "taqman_ref"


def run(sense, antisense, probe, show_structures=False):
    """
    print the values for the three oligos and their cross dimers
    when show_structures is set, also print an ASCII view of the winning self dimer, hairpin and cross dimer alignments under the table
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

    if show_structures:
        _print_structures(eng, seqs)

def _print_structures(eng, seqs, top_n=4):
    """
    ASCII rendering of the ranked alignments (best first)
    """
    print()
    print("structures")
    print("==========")
    for name, s in seqs.items():
        s = s.upper()
        print(f"\nself dimer  {name}")
        print(render.render_duplex_structures(s, s, eng.self_dimer_structures(s), top_n))
        print(f"\nhairpin  {name}")
        print(render.render_hairpin_structures(s, eng.hairpin_structures(s), top_n))
    for x, y in (("sense", "antisense"), ("sense", "probe"), ("antisense", "probe")):
        a, b = seqs[x].upper(), seqs[y].upper()
        print(f"\ncross dimer  {x} {y}")
        print(render.render_duplex_structures(a, b, eng.cross_dimer_structures(a, b), top_n))

def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    show_structures = False
    for flag in ("--show-structures", "-s"):
        if flag in argv:
            show_structures = True
            argv = [a for a in argv if a != flag]
    if len(argv) != 3:
        print("usage  python -m flare [--show-structures] SENSE ANTISENSE PROBE")
        return 1
    run(*argv, show_structures=show_structures)
    return 0


if __name__ == "__main__":
    sys.exit(main())
