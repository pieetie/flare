"""
ASCII rendering of secondary structures (self dimer, cross dimer, hairpin)
"""
from .params import COMPLEMENT

STRONG = '|'
WEAK = '¦'        # broken bar
BLANK = ' '
NOT_FOUND = 'Not Found'


def _is_wc(a, b):
    return COMPLEMENT.get(a) == b


def _fmt_dg(dg):
    return f"{round(dg, 1):.1f}"


def _footer(dg, width):
    return _fmt_dg(dg).rjust(width)


######
def render_duplex_structure(top, bottom, struct):
    top = top.upper()
    br = bottom.upper()[::-1]
    na, nb = len(top), len(br)
    offset = struct['offset']
    run = set(range(struct['start'], struct['start'] + struct['length']))
    top_pad = max(0, -offset)
    bot_pad = max(0, offset)
    width = max(top_pad + na, bot_pad + nb)

    syms = []
    for c in range(width):
        ti = c - top_pad
        bj = c - bot_pad
        if 0 <= ti < na and 0 <= bj < nb and _is_wc(top[ti], br[bj]):
            syms.append(STRONG if ti in run else WEAK)
        else:
            syms.append(BLANK)

    line_top = "5' " + " " * top_pad + top + " 3'"
    line_sym = "   " + "".join(syms)
    line_bot = "3' " + " " * bot_pad + br + " 5'"
    foot_w = max(len(line_top), len(line_bot))
    line_dg = _footer(struct['dg'], foot_w)
    return "\n".join(s.rstrip() for s in (line_top, line_sym, line_bot, line_dg))


def render_duplex_structures(top, bottom, structs, top_n=None):
    if not structs:
        return NOT_FOUND
    if top_n is not None:
        structs = structs[:top_n]
    return "\n\n".join(render_duplex_structure(top, bottom, s) for s in structs)


######

def render_hairpin_structure(seq, struct):
    seq = seq.upper()
    i, j, L = struct['i'], struct['j'], struct['L']
    split = (i + j + 1) // 2 # loop centre
    # an odd loop leaves a half-column offset -> nudge the 5' arm so the stem aligns
    top_lead = (i + j + 1) - 2 * split  # 0 or 1
    bot_lead = 0
    top_arm = seq[:split][::-1]
    bot_arm = seq[split:] 
    width = max(top_lead + len(top_arm), bot_lead + len(bot_arm))

    syms = []
    for c in range(width):
        tc, bc = c - top_lead, c - bot_lead
        ta = top_arm[tc] if 0 <= tc < len(top_arm) else None
        ba = bot_arm[bc] if 0 <= bc < len(bot_arm) else None
        if ta is not None and ba is not None and _is_wc(ta, ba):
            m_top = split - 1 - tc      # seq index of the 5' base in this column
            syms.append(STRONG if i <= m_top <= i + L - 1 else WEAK)
        else:
            syms.append(BLANK)

    line_top = "/" + " " * top_lead + top_arm + " 5'"
    line_sym = " " + "".join(syms)
    line_bot = "\\" + " " * bot_lead + bot_arm + " 3'"
    foot_w = max(len(line_top), len(line_sym), len(line_bot))
    line_dg = _footer(struct['dg'], foot_w)
    return "\n".join(s.rstrip() for s in (line_top, line_sym, line_bot, line_dg))

def render_hairpin_structures(seq, structs, top_n=None):
    if not structs:
        return NOT_FOUND
    if top_n is not None:
        structs = structs[:top_n]
    return "\n\n".join(render_hairpin_structure(seq, s) for s in structs)
