# RHUM-RUM periodic transient removal — version 2

Remove repeating instrument transients from a local MiniSEED stream using TiSKitPy.
The station periods and clipping levels are historical starting guesses, not
universal calibrations. Validate them against your instrument and waveform units.

## Run

Python 3.10 or newer:

```bash
git clone https://github.com/MohammadAmin-Aminian/Transients.git
cd Transients
python -m pip install -r requirements.txt
python Transients_Removal_24_Aug_23.py input.mseed cleaned.mseed \
  --station RR38 --transient-start 2012-10-12T00:00:00 --interactive
```

Replace the sample timestamp with the **observed transient onset**. Interactive
calibration opens plots requiring manual inspection and closure. Supported station
IDs are RR28, RR29, RR31, RR34, RR36, RR38, RR40, RR50 and RR52. Use `--help` for options.
Default windows are five days. Every chunk, including the final one, must contain
at least two transient periods; adjust `--window-days` if the tail is too short.

Input must contain one vertical channel and aligned, finite, gap-free traces with
identical start time, sample rate and sample count. Repair gaps and alignment
explicitly before running; the program does not silently interpolate missing data.
`--rotate` enables tilt correction and requires the horizontal components expected
by TiSKitPy (normally BH1/BH2). `--earthquake-spans` enables earthquake exclusion and
may contact a remote event service. Without these flags the workflow uses local data.
Existing output paths are refused. Output is float64 MiniSEED, retaining channel IDs.

## Fixes and validation

Version 2 uses the supported TiSKitPy 2.3.1 `PeriodicTransient` API, eliminates personal
absolute paths and unused station downloads, replaces the vertical trace in the
actual parent stream, appends each cleaned chunk once and preserves boundary samples.
No processing starts during import. Raw input files are never overwritten.

```bash
python -m pytest -q
```

Tests verify parent-trace replacement, timing validation and lossless chunk boundaries.
Scientific effectiveness on original survey waveforms remains unverified; inspect
before/after waveforms and spectra. Transient removal can also remove real signals
with similar timing, so calibration and earthquake exclusion matter.

Reference: [TiSKitPy periodic transients](https://tiskitpy.readthedocs.io/latest/periodic_transients.html).
Author: Mohammad Amin Aminian. No license was present in the original repository;
no additional reuse rights are asserted here.
