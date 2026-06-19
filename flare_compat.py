"""
Single-file oligo thermodynamics for a TaqMan assay, built on published nearest-neighbor science.

No dependencies beyond the Python standard library; download this file and call assay().
Parameters are calibrated against observed reference outputs for interoperability.

    from flare_compat import assay
    assay("TCTAACTAGCACACTAACTAATGTCA", "AACACTTGTGCGGTAACCTC", "CAGAATGTGTTAACCTGTCTTCT")

Defaults are the TaqMan reference conditions: 0.25 nM oligo, 50 mM monovalent, 5 mM free magnesium.
Parameters come from published literature; the calibrated values are tuned against observed reference outputs only.
See NOTICE.md for affiliation and trademark disclaimers.
"""
import math

R = 1.98720425864083
T37 = 310.15
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

# santalucia 1998 nearest neighbor dH (kcal/mol) and dS (cal/mol/K), used for tm
NN = {
    'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2),
    'AT': (-7.2, -20.4),
    'TA': (-7.2, -21.3),
    'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7),
    'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
    'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0),
    'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
    'CG': (-10.6, -27.2),
    'GC': (-9.8, -24.4),
    'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9),
}
INIT_GC = (0.1, -2.8)
INIT_AT = (2.3, 4.1)
SYMMETRY = (0.0, -1.4)

VON_AHSEN_K = 126.49
X_FACTOR = 1

# dimer dG37 per nearest neighbor step, calibrated fit, complement filled in
_DG37_RAW = {
    'AA': -1.1370, 'AT': -0.9730, 'TA': -0.7130, 'CA': -1.5870, 'GT': -1.5730,
    'CT': -1.3955, 'GA': -1.4430, 'CG': -2.3500, 'GC': -2.3830, 'GG': -1.9250,
}
DIMER_INIT = 1.87
DIMER_AT_PENALTY = 0.1455
DIMER_THRESHOLD = 0.35
DIMER_MIN_RUN = 3
HAIRPIN_NUCLEATION = 1.88
HAIRPIN_AT_PENALTY = 0.14
HAIRPIN_THRESHOLD = 0.0
HAIRPIN_MIN_STEM = 3
HAIRPIN_MIN_LOOP = 3


def _revcomp(seq):
    return ''.join(COMPLEMENT[b] for b in reversed(seq))


def _expand(base):
    out = dict(base)
    for k, v in base.items():
        out.setdefault(_revcomp(k), v)
    return out


DG37 = _expand(_DG37_RAW)


def na_total_mM(monovalent_mM, mg_free_mM):
    """
    total sodium equivalent in mM from monovalent and free magnesium
    """
    return monovalent_mM + VON_AHSEN_K * math.sqrt(mg_free_mM)


def gc_percent(seq):
    """
    percent of g and c bases
    """
    seq = seq.upper()
    return round(100.0 * (seq.count('G') + seq.count('C')) / len(seq), 2)


def gc_clamp(seq, window=5):
    """
    longest run of g or c inside the last few bases at the 3 prime end
    """
    seq = seq.upper()
    best = cur = 0
    for b in seq[-window:]:
        if b in 'GC':
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


def tm(seq, dna_nM=0.25, monovalent_mM=50.0, mg_free_mM=5.0):
    """
    melting temperature in celsius for the oligo against its complement
    """
    seq = seq.upper()
    na = na_total_mM(monovalent_mM, mg_free_mM) / 1000.0
    ct = dna_nM * 1e-9
    x = 1 if seq == _revcomp(seq) else X_FACTOR
    dH = dS = 0.0
    for i in range(len(seq) - 1):
        h, s = NN[seq[i:i + 2]]
        dH += h
        dS += s
    for end in (seq[0], seq[-1]):
        h, s = INIT_GC if end in 'GC' else INIT_AT
        dH += h
        dS += s
    if seq == _revcomp(seq):
        dH += SYMMETRY[0]
        dS += SYMMETRY[1]
    dS_salt = dS + 0.368 * (len(seq) - 1) * math.log(na)
    tm_k = (dH * 1000.0) / (dS_salt + R * math.log(ct / x))
    return round(tm_k - 273.15, 2)


def _report(dg, threshold):
    if dg is None or dg > threshold:
        return 0.0
    return round(dg, 1)


def _run_dg(top):
    if len(top) < 2:
        return None
    dg = DIMER_INIT + sum(DG37[top[i:i + 2]] for i in range(len(top) - 1))
    dg += DIMER_AT_PENALTY * ((top[0] in 'AT') + (top[-1] in 'AT'))
    return dg


def _best_duplex(a, b):
    a = a.upper()
    br = b.upper()[::-1]
    na, nb = len(a), len(br)
    best = None
    for d in range(-(nb - 1), na):
        run = []
        for j in range(nb):
            i = j + d
            if 0 <= i < na and COMPLEMENT.get(a[i]) == br[j]:
                run.append(a[i])
            else:
                if len(run) >= DIMER_MIN_RUN:
                    dg = _run_dg(''.join(run))
                    if dg is not None and (best is None or dg < best):
                        best = dg
                run = []
        if len(run) >= DIMER_MIN_RUN:
            dg = _run_dg(''.join(run))
            if dg is not None and (best is None or dg < best):
                best = dg
    return _report(best, DIMER_THRESHOLD)


def self_dimer(seq):
    """
    self dimer energy in kcal/mol
    """
    return _best_duplex(seq, seq)


def cross_dimer(a, b):
    """
    cross dimer energy in kcal/mol between two oligos
    """
    return _best_duplex(a, b)


def hairpin(seq):
    """
    hairpin energy in kcal/mol found by trying every loop then growing the stem
    """
    seq = seq.upper()
    n = len(seq)
    best = None
    for p in range(1, n):
        for q in range(p + HAIRPIN_MIN_LOOP - 1, n - 1):
            a, b = p - 1, q + 1
            stem5 = []
            while a >= 0 and b < n and COMPLEMENT.get(seq[a]) == seq[b]:
                stem5.append(seq[a])
                a -= 1
                b += 1
            if len(stem5) < HAIRPIN_MIN_STEM:
                continue
            top = ''.join(reversed(stem5))
            dg = HAIRPIN_NUCLEATION + sum(DG37[top[k:k + 2]] for k in range(len(top) - 1))
            dg += HAIRPIN_AT_PENALTY * ((top[0] in 'AT') + (top[-1] in 'AT'))
            if best is None or dg < best:
                best = dg
    return _report(best, HAIRPIN_THRESHOLD)


def analyze(seq, dna_nM=0.25, monovalent_mM=50.0, mg_free_mM=5.0):
    """
    full set of values for one oligo
    """
    seq = seq.upper()
    return {
        'sequence': seq,
        'length': len(seq),
        'tm_C': tm(seq, dna_nM, monovalent_mM, mg_free_mM),
        'gc_pct': gc_percent(seq),
        'gc_clamp': gc_clamp(seq),
        'self_dimer_dG_kcal': self_dimer(seq),
        'hairpin_dG_kcal': hairpin(seq),
    }


def assay(sense, antisense, probe, dna_nM=0.25, monovalent_mM=50.0, mg_free_mM=5.0):
    """
    full taqman assay, three oligos plus the three cross dimers
    """
    oligos = {
        'sense': analyze(sense, dna_nM, monovalent_mM, mg_free_mM),
        'antisense': analyze(antisense, dna_nM, monovalent_mM, mg_free_mM),
        'probe': analyze(probe, dna_nM, monovalent_mM, mg_free_mM),
    }
    cross = {
        'sense_antisense': cross_dimer(sense, antisense),
        'sense_probe': cross_dimer(sense, probe),
        'antisense_probe': cross_dimer(antisense, probe),
    }
    return {'oligos': oligos, 'cross_dimers': cross}


if __name__ == '__main__':
    import sys
    if len(sys.argv) != 4:
        print("usage  python flare_compat.py SENSE ANTISENSE PROBE")
        sys.exit(1)
    out = assay(*sys.argv[1:4])
    for name, r in out['oligos'].items():
        print(f"{name:10} Tm {r['tm_C']:6.2f}  GC% {r['gc_pct']:5.1f}  clamp {r['gc_clamp']}"
              f"  self {r['self_dimer_dG_kcal']:+.1f}  hairpin {r['hairpin_dG_kcal']:+.1f}")
    for k, v in out['cross_dimers'].items():
        print(f"cross {k} {v:+.1f}")
