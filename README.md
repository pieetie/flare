# Flare

Flare is an independent oligo thermodynamics engine for TaqMan assays, built on published nearest-neighbor science.

Parameters come from published science. The optional `calibrated` mode tunes a handful of parameters against observed reference values, with no access to any third-party source code or binaries. See `NOTICE.md`.

![Flare error per metric against recorded reference values](figures/flare_validation.png)

## Use

Run Flare on an assay:

    uv run python -m flare SENSE ANTISENSE PROBE

Save the comparison figure as PNG:

    uv run python -m flare.plots

Run the frozen tests:

    uv run python -m pytest tests/ -q
