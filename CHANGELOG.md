# Changelog

## v1.1.0 — 2026-06-25

- ASCII visualization of secondary structures (self dimer, cross dimer, hairpin)
- New `--show-structures` / `-s` CLI flag to print ranked alignments under the table
- Engine helpers `self_dimer_structures`, `cross_dimer_structures` and `hairpin_structures` for ranked alignments (visualization only, no effect on dG)

## v1.0.0 — 2026-06-19

- Release of Flare thermodynamics engine for TaqMan assays
- Simple CLI, validation plots and pytest suite with golden master
- Published nearest-neighbor parameters; optional calibrated mode
