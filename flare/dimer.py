"""
self dimer, cross dimer and hairpin energies by nearest neighbor stacking
"""
from .params import NN_SETS, COMPLEMENT

T37 = 310.15


def dG37_table(nn_set, temp_K=T37):
    """
    per step dG37 in kcal/mol from dH and dS
    """
    return {k: h - temp_K * s / 1000.0 for k, (h, s) in nn_set['nn'].items()}


def get_dG37(config, temp_K=T37):
    """
    use the fitted dG37 table if given else build one from the nn set
    """
    direct = config.get('dimer_dG37_table')
    if direct is not None:
        return direct
    return dG37_table(NN_SETS[config['nn_params']], temp_K)


def is_wc(a, b):
    return COMPLEMENT.get(a) == b


def _run_dG(top_run, dG_nn, init, at_penalty=0.0):
    """
    energy of one contiguous paired run given its top strand
    """
    if len(top_run) < 2:
        return None
    dg = init
    for i in range(len(top_run) - 1):
        dg += dG_nn[top_run[i:i + 2]]
    if at_penalty:
        dg += at_penalty * ((top_run[0] in 'AT') + (top_run[-1] in 'AT'))
    return dg


def best_duplex(a, b, config, temp_K=T37):
    """
    most stable antiparallel duplex between two strands given 5 to 3
    """
    dG_nn = get_dG37(config, temp_K)
    init = config.get('dimer_init_dG', 1.96)
    at_pen = config.get('dimer_terminal_at_penalty', 0.0)
    min_run = config.get('dimer_min_run', 2)
    a = a.upper()
    br = b.upper()[::-1]
    na, nb = len(a), len(br)
    best = (None, None)

    def consider(run, d, col_start):
        if len(run) < min_run:
            return None
        dg = _run_dG(''.join(run), dG_nn, init, at_pen)
        if dg is not None and (best[0] is None or dg < best[0]):
            return dg, dict(offset=d, top=''.join(run), start=col_start)
        return None

    for d in range(-(nb - 1), na):
        run = []
        col_start = None
        for j in range(nb):
            i = j + d
            if 0 <= i < na and is_wc(a[i], br[j]):
                if not run:
                    col_start = i
                run.append(a[i])
            else:
                r = consider(run, d, col_start)
                if r:
                    best = r
                run = []
        r = consider(run, d, col_start)
        if r:
            best = r
    return (best[0] if best[0] is not None else 0.0, best[1])


def self_dimer(seq, config):
    return best_duplex(seq, seq, config)


def cross_dimer(a, b, config):
    return best_duplex(a, b, config)


def hairpin(seq, config, min_loop=3, temp_K=T37):
    """
    most stable hairpin found by trying every loop then growing the stem outward
    """
    dG_nn = get_dG37(config, temp_K)
    nucleation = config.get('hairpin_nucleation_dG', 1.8)
    at_pen = config.get('hairpin_terminal_at_penalty', 0.0)
    min_stem = config.get('hairpin_min_stem', 2)
    seq = seq.upper()
    n = len(seq)
    best = (0.0, None)
    for p in range(1, n):
        for q in range(p + min_loop - 1, n - 1):
            a, b = p - 1, q + 1
            stem5 = []
            while a >= 0 and b < n and is_wc(seq[a], seq[b]):
                stem5.append(seq[a])
                a -= 1
                b += 1
            L = len(stem5)
            if L < min_stem:
                continue
            top = ''.join(reversed(stem5))
            dg = nucleation
            for k in range(L - 1):
                dg += dG_nn[top[k:k + 2]]
            if at_pen:
                dg += at_pen * ((top[0] in 'AT') + (top[-1] in 'AT'))
            if dg < best[0]:
                best = (dg, dict(i=a + 1, j=b - 1, L=L, loop=q - p + 1, stem=top))
    return best


def report_dG(dg, threshold=-0.1, ndp=1):
    """
    round to one decimal and report 0.0 when weaker than the threshold
    """
    if dg is None or dg > threshold:
        return 0.0
    return round(dg, ndp)
