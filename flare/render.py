"""
ASCII rendering of secondary structures (self dimer, cross dimer, hairpin)
"""
from .params import COMPLEMENT

STRONG = '|'
WEAK = '¦' # broken bar, the "weak" pairing glyph
BLANK = ' '
STRONG_PAIRS = {('G', 'C'), ('C', 'G')}
NOT_FOUND = 'Not Found'

def _is_wc(a, b):
    return COMPLEMENT.get(a) == b

def _pair_symbol(a, b):
    if not _is_wc(a, b):
        return BLANK
    return STRONG if (a, b) in STRONG_PAIRS else WEAK

def _fmt_dg(dg):
    return f"{dg:.1f}"

def _footer(dg, width):
    return _fmt_dg(dg).rjust(width)

def render_duplex(top, bottom, offset, dg):
    """
    Three line antiparallel duplex render plus a right aligned dG footer.

    top      : the 5'->3' strand drawn on the first line
    bottom   : the partner strand given 5'->3'; it is drawn 3'->5' (reversed)
               so paired bases line up vertically
    offset   : the winning relative shift the engine found. column i of the top
               pairs column (i - offset) of the reversed bottom. positive offset
               shifts the bottom strand to the right, negative to the left.
    dg       : the scalar dG already reported by the engine (kcal/mol)
    """
    top = top.upper()
    br = bottom.upper()[::-1] # bottom written 3'->5'
    na, nb = len(top), len(br)
    top_pad = max(0, -offset)
    bot_pad = max(0, offset)
    width = max(top_pad + na, bot_pad + nb) 
    syms = []
    for c in range(width):
        ti = c - top_pad
        bj = c - bot_pad
        if 0 <= ti < na and 0 <= bj < nb: 
            syms.append(_pair_symbol(top[ti], br[bj]))
        else:
            syms.append(BLANK)

    line_top = "5' " + " " * top_pad + top + " 3'" 
    line_sym = "   " + "".join(syms)
    line_bot = "3' " + " " * bot_pad + br + " 5'"
    foot_w = max(len(line_top), len(line_bot))
    line_dg = _footer(dg, foot_w)  

    return "\n".join(s.rstrip() for s in (line_top, line_sym, line_bot, line_dg)) 

def render_hairpin(seq, st, dg):
    seq = seq.upper()
    i, j, L = st['i'], st['j'], st['L']
    stem5 = seq[i:i + L] 
    arm3 = seq[j - L + 1:j + 1] 
    arm3_disp = arm3[::-1]
    loop = seq[i + L:j - L + 1] 

    syms = "".join(_pair_symbol(stem5[k], arm3_disp[k]) for k in range(L))

    line_top = "5' " + stem5 + " \\" 
    line_sym = "   " + syms + "  ) " + loop
    line_bot = "3' " + arm3_disp + " /" 
    foot_w = max(len(line_top), len(line_sym), len(line_bot))
    line_dg = _footer(dg, foot_w) 

    return "\n".join(s.rstrip() for s in (line_top, line_sym, line_bot, line_dg))

def render_self_dimer(seq, dg, st):
    """
    Self dimer: both strands are the same sequence
    'Not Found' if no structure
    """
    if st is None or dg == 0.0:
        return NOT_FOUND
    return render_duplex(seq, seq, st['offset'], dg) 
  
def render_cross_dimer(a, b, dg, st):
    """
    Cross dimer between strands a and b
    'Not Found' if no structure
    """
    if st is None or dg == 0.0:
        return NOT_FOUND
    return render_duplex(a, b, st['offset'], dg)

def render_hairpin_structure(seq, dg, st):
    """
    Hairpin of one strand
    'Not Found' if no structure
    """
    if st is None or dg == 0.0:
        return NOT_FOUND
    return render_hairpin(seq, st, dg) 