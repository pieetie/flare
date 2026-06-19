"""
golden master, compares every flare value to the frozen snapshot
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from flare import loader                       # noqa: E402
from flare.core import Engine                  # noqa: E402
import tests.gen_golden as gen                  # noqa: E402

GOLDEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'golden_engine.json')


def test_golden_master_unchanged():
    with open(GOLDEN) as f:
        frozen = json.load(f)
    current = gen.build()
    diffs = []

    def walk(path, a, b):
        if isinstance(a, dict):
            if not isinstance(b, dict) or set(a) != set(b):
                diffs.append(f"{path} dict keys differ")
                return
            for k in a:
                walk(f"{path}.{k}", a[k], b[k])
        elif isinstance(a, list):
            if not isinstance(b, list) or len(a) != len(b):
                diffs.append(f"{path} list length changed")
                return
            for i, (x, y) in enumerate(zip(a, b)):
                walk(f"{path}[{i}]", x, y)
        else:
            if a != b:
                diffs.append(f"{path} {b!r} now vs {a!r} frozen")

    walk("root", frozen, current)
    assert not diffs, "golden master drift (%d)\n  %s" % (len(diffs), "\n  ".join(diffs[:40]))


def test_golden_file_is_complete():
    # catch a truncated or empty snapshot
    with open(GOLDEN) as f:
        frozen = json.load(f)
    assert len(frozen['oligos']) == 93, len(frozen['oligos'])
    assert len(frozen['cross_dimers']) == 93, len(frozen['cross_dimers'])
    assert len(frozen['arbitrary']) == len(gen.ARBITRARY)
    assert len(frozen['arbitrary_cross']) == len(gen.ARBITRARY) * (len(gen.ARBITRARY) - 1) // 2


if __name__ == '__main__':
    for name, fn in sorted(globals().items()):
        if name.startswith('test_'):
            try:
                fn()
                print("PASS ", name)
            except AssertionError as e:
                print("FAIL ", name, ":", str(e)[:500])
