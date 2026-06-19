"""
default model settings used by the engine
"""

DEFAULTS = {
    "nn_params": "santalucia_2004",
    "salt_correction": "santalucia_1998",
    "von_ahsen_k": 126.49,
    "x_factor": 1,
    "dimer_init_dG": 1.65,
    "dimer_threshold": -0.65,
    "hairpin_nucleation_dG": 1.8,
}


def load_config():
    """
    return a fresh copy of the default settings
    """
    return dict(DEFAULTS)
