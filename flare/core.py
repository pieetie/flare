"""
Engine that turns an oligo and reaction conditions into thermodynamic values.
"""
from . import loader, salt as saltmod, tm as tmmod, dimer, params
from .config import load_config


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


class Engine:
    """
    Compute Tm, GC and dimer energies in literature or exact compatibility mode.
    """

    def __init__(self, config=None, calibration_mode='literature', reference_cache=False):
        self.cfg = dict(config or load_config())
        self.mode = calibration_mode
        self.conditions = loader.load_conditions()
        self.reference = None
        if calibration_mode == 'beacon_exact':
            self.cfg['dimer_dG37_table'] = params.DG37_BEACON_FIT
            self.cfg['dimer_init_dG'] = params.DG37_BEACON_FIT_INIT
            self.cfg['dimer_terminal_at_penalty'] = params.DG37_BEACON_FIT_AT_PENALTY
            self.cfg['hairpin_nucleation_dG'] = params.DG37_BEACON_FIT_HAIRPIN_NUCLEATION
            self.cfg['hairpin_terminal_at_penalty'] = params.DG37_BEACON_FIT_HAIRPIN_AT_PENALTY
            self.cfg['dimer_min_run'] = 3
            self.cfg['hairpin_min_stem'] = 3
            self.cfg['dimer_threshold'] = params.DG37_BEACON_FIT_THRESHOLD
            self.cfg['hairpin_threshold'] = params.DG37_BEACON_FIT_HAIRPIN_THRESHOLD
        elif calibration_mode != 'literature':
            raise ValueError(f"unknown calibration_mode {calibration_mode!r}")
        if reference_cache:
            from .reference import Reference
            self.reference = Reference()

    def _cond(self, conditions):
        if isinstance(conditions, str):
            return self.conditions[conditions]
        return conditions

    def na_total_mM(self, conditions):
        c = self._cond(conditions)
        k = self.cfg.get('von_ahsen_k', saltmod.VON_AHSEN_K)
        return saltmod.mg_to_na_total_mM(c['monovalent_mM'], c['mg_free_mM'], k)

    def tm(self, seq, conditions):
        c = self._cond(conditions)
        if self.reference is not None and isinstance(conditions, str):
            hit = self.reference.tm(seq, conditions)
            if hit is not None:
                return hit
        ct = c['nucleic_acid_nM'] * 1e-9
        na = self.na_total_mM(c) / 1000.0
        return round(tmmod.tm(seq, ct, na, self.cfg), 2)

    def self_dimer_dG(self, seq):
        if self.reference is not None:
            hit = self.reference.self_dimer(seq)
            if hit is not None:
                return hit, None
        dg, st = dimer.self_dimer(seq, self.cfg)
        return dimer.report_dG(dg, self.cfg.get('dimer_threshold', -0.45)), st

    def cross_dimer_dG(self, a, b):
        if self.reference is not None:
            hit = self.reference.cross_dimer(a, b)
            if hit is not None:
                return hit, None
        dg, st = dimer.cross_dimer(a, b, self.cfg)
        return dimer.report_dG(dg, self.cfg.get('dimer_threshold', -0.45)), st

    def hairpin_dG(self, seq):
        if self.reference is not None:
            hit = self.reference.hairpin(seq)
            if hit is not None:
                return hit, None
        dg, st = dimer.hairpin(seq, self.cfg)
        thr = self.cfg.get('hairpin_threshold', self.cfg.get('dimer_threshold', -0.45))
        return dimer.report_dG(dg, thr), st

    def analyze_oligo(self, seq, conditions):
        """
        full set of values for one oligo
        """
        seq = seq.upper()
        sd, _ = self.self_dimer_dG(seq)
        hp, _ = self.hairpin_dG(seq)
        return dict(
            sequence=seq, length=len(seq),
            tm_C=self.tm(seq, conditions),
            gc_pct=gc_percent(seq), gc_clamp=gc_clamp(seq),
            self_dimer_dG_kcal=sd, hairpin_dG_kcal=hp,
        )
