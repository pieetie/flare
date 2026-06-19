"""
salt handling, magnesium to sodium equivalent and the tm salt term
"""
import math

VON_AHSEN_K = 126.49


def mg_to_na_total_mM(monovalent_mM, mg_free_mM, k=VON_AHSEN_K):
    """
    total sodium equivalent in mM from monovalent and free magnesium
    """
    return monovalent_mM + k * math.sqrt(mg_free_mM)


def salt_santalucia_1998_dS(dS, na_molar, n_phosphates):
    """
    santalucia 1998 entropy correction for sodium
    """
    return dS + 0.368 * (n_phosphates - 1) * math.log(na_molar)
