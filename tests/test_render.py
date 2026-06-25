"""
frozen master for the ASCII structure renderer
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from flare.core import Engine          # noqa: E402
from flare import render               # noqa: E402

# the broken bar weak-pair glyph, spelled out so the file stays ascii-safe
W = "¦"


def _eng():
    return Engine(calibration_mode="calibrated")


def test_self_dimer_example_minus_0_9():
    """the user example: self dimer of the sense oligo, dG -0.9, offset 0."""
    e = _eng()
    seq = "CAGTCAGCGGGCAATGAATC"
    dg, st = e.self_dimer_dG(seq)
    assert dg == -0.9
    assert st["offset"] == 0
    out = render.render_self_dimer(seq, dg, st)
    expected = "\n".join([
        "5' CAGTCAGCGGGCAATGAATC 3'",
        "    " + W + " " + W + "|" + W + "  |  |  " + W + "|" + W + " " + W,
        "3' CTAAGTAACGGGCGACTGAC 5'",
        "                      -0.9",
    ])
    assert out == expected, repr(out)


def test_self_dimer_offset_zero_full_match():
    """offset 0, every column pairs; sanity on perfect antiparallel duplex."""
    e = _eng()
    seq = "ACGTACGT"
    dg, st = e.self_dimer_dG(seq)
    assert st["offset"] == 0
    out = render.render_self_dimer(seq, dg, st)
    lines = out.splitlines()
    assert lines[0] == "5' ACGTACGT 3'"
    assert lines[2] == "3' TGCATGCA 5'"
    # 8 columns: A:T weak, C:G/G:C strong, by chemistry lookup
    assert lines[1] == "   " + W + "||" + W + W + "||" + W
    # footer carries the engine dG, right aligned to the strand width
    assert lines[3].strip() == f"{dg:.1f}"
    assert lines[3].endswith(f"{dg:.1f}")
    assert len(lines[3]) == len(lines[0])


def test_cross_dimer_negative_offset():
    """negative offset must pad the top strand, not the bottom."""
    e = _eng()
    a, b = "GGGGGGGG", "CCCCAAAA"
    dg, st = e.cross_dimer_dG(a, b)
    assert st["offset"] == -4
    out = render.render_cross_dimer(a, b, dg, st)
    expected = "\n".join([
        "5'     GGGGGGGG 3'",
        "       ||||",
        "3' AAAACCCC 5'",
        "              " + f"{dg:.1f}",
    ])
    assert out == expected, repr(out)
    # the four bonds sit exactly above the paired C's of the bottom strand
    top, sym, bot, _ = out.splitlines()
    first_bond = sym.index("|")
    assert top[first_bond] == "G"
    assert bot[first_bond] == "C"


def test_hairpin_render():
    """hairpin draws two arms of one strand plus the loop bases."""
    e = _eng()
    seq = "CAGTCAGCGGGCAATGAATC"
    dg, st = e.hairpin_dG(seq)
    assert dg == -0.9
    out = render.render_hairpin_structure(seq, dg, st)
    expected = "\n".join([
        "5' TCA \\",
        "   " + W + "|" + W + "  ) GCGGGCAA",
        "3' AGT /",
        "              -0.9",
    ])
    assert out == expected, repr(out)


def test_not_found_when_no_structure():
    """no pairing under threshold renders 'Not Found', never an empty diagram."""
    e = _eng()
    # short A-only oligo: engine returns dg 0.0 and st None
    dg, st = e.self_dimer_dG("AAAAA")
    assert dg == 0.0 and st is None
    assert render.render_self_dimer("AAAAA", dg, st) == "Not Found"
    # cross dimer with no complementarity
    dg, st = e.cross_dimer_dG("AAAA", "AAAA")
    assert render.render_cross_dimer("AAAA", "AAAA", dg, st) == "Not Found"
    # a probe that is below threshold both ways
    dg, st = e.hairpin_dG("AAAAA")
    assert render.render_hairpin_structure("AAAAA", dg, st) == "Not Found"


def test_columns_align_under_top_strand():
    """every symbol must sit on a column where both strands carry a base."""
    e = _eng()
    for seq in ("CAGTCAGCGGGCAATGAATC", "ACGTACGT", "GCATGCATGCAT"):
        dg, st = e.self_dimer_dG(seq)
        if st is None:
            continue
        out = render.render_self_dimer(seq, dg, st)
        top, sym, bot, _ = out.splitlines()
        for c, ch in enumerate(sym):
            if ch in ("|", W):
                assert c < len(top) and top[c] not in (" ", "'")
                assert c < len(bot) and bot[c] not in (" ", "'")


def test_renderer_does_not_change_engine_dg():
    """the footer dG is exactly the engine scalar, formatted, nothing recomputed."""
    e = _eng()
    seq = "CAGTCAGCGGGCAATGAATC"
    dg, st = e.self_dimer_dG(seq)
    out = render.render_self_dimer(seq, dg, st)
    assert out.splitlines()[-1].strip() == f"{dg:.1f}"


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
