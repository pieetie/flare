import os
import sys

# Ensure the repo root is importable (flare, tests) under any pytest invocation.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
