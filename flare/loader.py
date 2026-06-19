"""
read oligos, conditions and assays from the data folder
"""
import os
import yaml

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')


def _load(rel):
    with open(os.path.join(DATA, rel)) as f:
        return yaml.safe_load(f)


def load_index():
    return _load('index.yaml')


def load_oligos():
    """
    map oligo id to sequence
    """
    out = {}
    for entry in load_index()['oligos']:
        d = _load(entry['file'])
        out[d['id']] = d['sequence'].upper()
    return out


def load_conditions():
    """
    map condition id to concentrations
    """
    out = {}
    for entry in load_index()['conditions']:
        d = _load(entry['file'])
        out[d['id']] = d['concentrations']
    return out


def load_assays():
    """
    list of assay records
    """
    return [_load(e['file']) for e in load_index()['assays']]


def load_structures(rel):
    return _load(rel)
