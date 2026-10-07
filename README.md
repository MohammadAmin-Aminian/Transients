# OBS Transient Cleaner

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

**Periodic instrument-transient detection and removal for ocean-bottom seismic data.**

[![Regression tests](https://github.com/MohammadAmin-Aminian/obs-transient-cleaner/actions/workflows/tests.yml/badge.svg)](https://github.com/MohammadAmin-Aminian/obs-transient-cleaner/actions/workflows/tests.yml)
[![TiSKitPy](https://img.shields.io/badge/TiSKitPy-2.3.1-blue)](https://tiskitpy.readthedocs.io/)

This repository provides a focused workflow for removing repeating instrumental transients from ocean-bottom seismometer records. It was developed around the RHUM-RUM experiment, where approximately hourly mass-positioning events can contaminate long-period seismic processing.

The code is intentionally narrow in scope: it handles **periodic transient suppression**, timing integrity, chunked processing and optional tilt-related rotation before cleaning. It is not a general denoising package.

## Scientific motivation

Long-period OBS processing is especially sensitive to non-seismic signals. In broadband seafloor records, repeating instrument events can leak into spectral estimates, coherence calculations and later compliance processing if they are not identified and removed consistently.

The RHUM-RUM compliance workflow used careful preprocessing before estimating the relationship between vertical displacement and pressure. This repository isolates one of those preprocessing tasks so it can be tested and reused independently.

Associated scientific context: [Aminian et al. (2025), Geophysical Journal International](https://doi.org/10.1093/gji/ggaf253).

## What the tool does

```text
MiniSEED input
    |
    v
validate timing / channels / gaps
    |
    v
split into processing windows
    |
    +--> optional tilt-related rotation
    |
    v
construct periodic transient model
    |
    v
estimate transient waveform
    |
    v
remove repeating transient
    |
    v
merge chunks without losing boundaries
    |
    v
cleaned MiniSEED output
```

Key features:

- station-specific historical transient periods and clipping ranges;
- explicit transient onset supplied by the user;
- deterministic chunking with no duplicated or dropped boundary samples;
- preservation of trace identity and timing;
- optional `CleanRotator` step through TiSKitPy;
- optional earthquake-span exclusion during calibration;
- synthetic integration test using the real TiSKitPy transient remover;
- refusal to overwrite existing outputs.

## Supported RHUM-RUM stations

Current historical presets are included for:

`RR28`, `RR29`, `RR31`, `RR34`, `RR36`, `RR38`, `RR40`, `RR50`, `RR52`.

These values are **starting parameters from the original workflow**, not universal instrument calibrations. They should be inspected and recalibrated for other deployments or instruments.

## Installation

Python 3.10 or newer:

```bash
git clone https://github.com/MohammadAmin-Aminian/obs-transient-cleaner.git
cd obs-transient-cleaner
python -m pip install -e '.[dev]'
```

Dependencies include NumPy, ObsPy and TiSKitPy 2.3.1.

## Basic use

```bash
obs-transient-clean input.mseed cleaned.mseed \
    --station RR38 \
    --transient-start 2012-10-12T00:00:00
```

For interactive timing refinement:

```bash
obs-transient-clean input.mseed cleaned.mseed \
    --station RR38 \
    --transient-start 2012-10-12T00:00:00 \
    --interactive
```

Optional processing controls:

- `--window-days` — chunk duration, default 5 days;
- `--rotate` — apply TiSKitPy `CleanRotator` before transient removal;
- `--earthquake-spans` — enable earthquake exclusion during calibration;
- `--interactive` — inspect/refine transient timing.

The historical script entry point remains available as `python Transients_Removal_24_Aug_23.py ...` for backward compatibility. Run `obs-transient-clean --help` for the full command-line interface.

## Input requirements

The input stream must:

- contain exactly one vertical channel;
- contain traces from the requested station only;
- have aligned start times;
- have identical sampling rates and sample counts;
- be finite and gap-free.

The program deliberately rejects ambiguous or damaged inputs instead of silently interpolating them.

Each processing chunk must contain at least two transient periods. If the final chunk is too short, increase or decrease `--window-days`.

## Validation

Run:

```bash
python -m pytest -q
```

The current tests verify:

- exact chunk-boundary preservation;
- parent-stream trace replacement;
- rejection of misaligned traces;
- one-pass processing of each chunk;
- timing and identity preservation;
- substantial suppression of known synthetic periodic pulses using real TiSKitPy processing.

The synthetic benchmark requires the residual pulse RMS error to fall below 20% of the original pulse RMS for the tested configuration.

## Interpretation and limitations

Transient removal can also suppress real signals if they overlap the learned transient timing or waveform. For that reason:

- inspect before/after waveforms;
- compare power spectra before and after cleaning;
- exclude earthquakes when appropriate;
- do not transfer station-specific timing parameters blindly;
- validate the cleaned trace before downstream compliance or noise analysis.

This project does **not** claim to distinguish all instrumental artefacts from geophysical signals. It targets a specific class of repeating OBS transients.

## Relationship to ComPy

This repository is a standalone preprocessing tool.

- **OBS Transient Cleaner**: repeating instrumental transient removal.
- **ComPy**: broader compliance processing, calibration and inversion.

The functionality was used in the wider scientific workflow, but it is useful independently for OBS quality control and long-period preprocessing.

## Reference

Aminian, M. A., Crawford, W., Stutzmann, É., Montagner, J.-P., Cannat, M., & Hadziioannou, C. (2025). *Shallow crustal structures of the Indian ocean derived from compliance function analysis*. Geophysical Journal International, 242(3), ggaf253. https://doi.org/10.1093/gji/ggaf253

TiSKitPy periodic transient documentation: https://tiskitpy.readthedocs.io/latest/periodic_transients.html

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md) and [CHANGELOG.md](CHANGELOG.md).

**Author:** Mohammad Amin Aminian

## Related research software

This repository is part of a broader seismic/geophysical software portfolio:

- [ComPy](https://github.com/MohammadAmin-Aminian/ComPy) — seafloor compliance processing, DPG calibration and layered elastic inversion.
- [OBS Transient Cleaner](https://github.com/MohammadAmin-Aminian/obs-transient-cleaner) — periodic OBS instrument-transient removal.
- [ComPy Inversion Tuner](https://github.com/MohammadAmin-Aminian/compy-inversion-tuner) — reproducible tuning of compliance-inversion controls.
- [RHUM-RUM Geospatial Mapper](https://github.com/MohammadAmin-Aminian/rhum-rum-geospatial-mapper) — bathymetry, OBS-network and tectonic-context mapping.
- [VRE Seismic Enhancement](https://github.com/MohammadAmin-Aminian/vre-seismic-enhancement) — Virtual Resolution Enhancement for seismic sections.
- [Gabor Seismic Filter](https://github.com/MohammadAmin-Aminian/gabor-seismic-filter) — orientation-selective 2-D seismic filtering in MATLAB.


## License

This software is released under the **GNU General Public License v3.0 only (GPL-3.0-only)**. See [LICENSE](LICENSE) for the full terms.
