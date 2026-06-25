"""
frozen master for the ASCII structure renderer
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from flare.core import Engine          # noqa: E402
from flare import render               # noqa: E402

W = "¦"   # the broken-bar weak/incidental glyph


def _eng():
    return Engine(calibration_mode="calibrated")


def test_self_dimer_run_incident_shading():
    e = _eng()
    seq = "CAGTCAGCGGGCAATGAATC"
    sts = e.self_dimer_structures(seq)
    out = render.render_duplex_structure(seq, seq, sts[0])
    expected = "\n".join([
        "5' CAGTCAGCGGGCAATGAATC 3'",
        "    " + W + " |||  " + W + "  " + W + "  " + W + W + W + " " + W,
        "3' CTAAGTAACGGGCGACTGAC 5'",
        "                      -0.9",
    ])
    assert out == expected, repr(out)
    assert sts[0]["length"] == 3
    assert seq[sts[0]["start"]:sts[0]["start"] + 3] == "TCA"


def test_cross_dimer_matches_beacon_primary():
    e = _eng()
    a, b = "CAGTCAGCGGGCAATGAATC", "TCATCCTCAGTCCTTGCTTCT"
    sts = e.cross_dimer_structures(a, b)
    out = render.render_duplex_structure(a, b, sts[0])
    expected = "\n".join([
        "5' CAGTCAGCGGGCAATGAATC 3'",
        "             |||| " + W + W + " " + W,
        "3'       TCTTCGTTCCTGACTCCTACT 5'",
        "                             -3.1",
    ])
    assert out == expected, repr(out)


def test_cross_dimer_negative_offset_pads_top():
    e = _eng()
    a, b = "GGGGGGGG", "CCCCAAAA"
    sts = e.cross_dimer_structures(a, b)
    assert sts[0]["offset"] == -4
    out = render.render_duplex_structure(a, b, sts[0])
    top, sym, bot, foot = out.splitlines()
    assert top == "5'     GGGGGGGG 3'"
    assert bot == "3' AAAACCCC 5'"
    assert sym == "       ||||"
    assert foot.strip() == "-3.9"
    first = sym.index("|")
    assert top[first] == "G" and bot[first] == "C"


def test_hairpin_reference_style():
    e = _eng()
    seq = "CAGTCAGCGGGCAATGAATC"
    sts = e.hairpin_structures(seq)
    out = render.render_hairpin_structure(seq, sts[0])
    expected = "\n".join([
        "/GGCGACTGAC 5'",
        "  " + W + "  ||| " + W,
        "\\GCAATGAATC 3'",
        "          -0.9",
    ])
    assert out == expected, repr(out)


def test_hairpin_odd_loop_keeps_stem_aligned():
    e = _eng()
    seq = "TTCCACCCTTTCCTTCTGGT"
    for st in e.hairpin_structures(seq):
        out = render.render_hairpin_structure(seq, st)
        bars = out.splitlines()[1].count("|")
        assert bars == st["L"], (bars, st)


def test_not_found_when_no_structure():
    e = _eng()
    assert e.self_dimer_structures("TCATCCTCAGTCCTTGCTTCT") == []
    assert render.render_duplex_structures("TCATCCTCAGTCCTTGCTTCT",
                                           "TCATCCTCAGTCCTTGCTTCT", []) == "Not Found"
    assert render.render_hairpin_structures("AAAAA", []) == "Not Found"


def test_multi_structure_block_is_best_first():
    e = _eng()
    a, b = "CAGTCAGCGGGCAATGAATC", "TCATCCTCAGTCCTTGCTTCT"
    sts = e.cross_dimer_structures(a, b)
    block = render.render_duplex_structures(a, b, sts, top_n=4)
    dgs = [float(line) for line in block.splitlines()
           if line.strip().lstrip("-").replace(".", "").isdigit()]
    assert dgs == sorted(dgs)            # ascending = best first
    assert dgs[0] == -3.1


def test_footer_is_engine_dg_not_recomputed():
    e = _eng()
    seq = "CAGTCAGCGGGCAATGAATC"
    st = e.self_dimer_structures(seq)[0]
    out = render.render_duplex_structure(seq, seq, st)
    assert out.splitlines()[-1].strip() == f"{round(st['dg'], 1):.1f}"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    bad = 0
    for fn in fns:
        try:
            fn()
            print("PASS ", fn.__name__)
        except AssertionError as ex:
            bad += 1
            print("FAIL ", fn.__name__, ":", ex)
    print(f"\n{len(fns) - bad}/{len(fns)} passed")
    sys.exit(1 if bad else 0)
