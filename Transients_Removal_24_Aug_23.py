"""Remove calibrated periodic transients from local RHUM-RUM MiniSEED files."""

import argparse
from pathlib import Path

import numpy as np
from obspy import Stream, UTCDateTime, read

PARAMETERS = {
    "RR28": (3620.3, (-820, 330)),
    "RR29": (3620.3, (-290, 135)),
    "RR31": (3620.3, (-290, 135)),
    "RR34": (3620.3, (-290, 135)),
    "RR36": (3620.3, (-290, 135)),
    "RR38": (3620, (-235, 180)),
    "RR40": (3620.3, (-290, 135)),
    "RR50": (3620.3, (-290, 135)),
    "RR52": (3619.7, (-290, 135)),
}


def split_stream(stream, seconds):
    """Partition aligned traces with no duplicated or discarded boundary samples."""
    if not stream or not np.isfinite(seconds) or seconds <= 0:
        raise ValueError("Need a nonempty stream and positive window duration")
    reference = stream[0].stats
    for trace in stream:
        if (
            trace.stats.starttime != reference.starttime
            or trace.stats.sampling_rate != reference.sampling_rate
            or trace.stats.npts != reference.npts
        ):
            raise ValueError(
                "Traces must have identical start times, rates and lengths"
            )
        if np.ma.isMaskedArray(trace.data) or not np.isfinite(trace.data).all():
            raise ValueError("Gapped or nonfinite traces must be repaired explicitly")
    count = int(seconds * reference.sampling_rate)
    if count < 1 or reference.npts < 2:
        raise ValueError("Window or trace is too short")
    for offset in range(0, reference.npts, count):
        chunk = Stream()
        for trace in stream:
            part = trace.copy()
            part.data = trace.data[offset : offset + count].copy()
            part.stats.starttime = (
                trace.stats.starttime + offset / trace.stats.sampling_rate
            )
            chunk.append(part)
        yield chunk


def clean_chunk(chunk, cleaner):
    """Replace the vertical trace in the parent stream, keeping other channels."""
    result = chunk.copy()
    indices = [i for i, t in enumerate(result) if t.stats.channel.endswith("Z")]
    if len(indices) != 1:
        raise ValueError("Expected exactly one vertical channel")
    index = indices[0]
    original = result[index]
    cleaned = cleaner(original.copy())
    if (
        cleaned.stats.npts != original.stats.npts
        or cleaned.stats.starttime != original.stats.starttime
        or cleaned.stats.sampling_rate != original.stats.sampling_rate
    ):
        raise ValueError("Cleaner changed trace timing or length")
    cleaned.stats.network = original.stats.network
    cleaned.stats.station = original.stats.station
    cleaned.stats.location = original.stats.location
    cleaned.stats.channel = original.stats.channel
    result[index] = cleaned
    return result


def process(
    stream,
    station,
    transient_start,
    *,
    window_seconds=432000,
    rotate=False,
    interactive=False,
    earthquake_spans=False,
):
    if station not in PARAMETERS:
        raise ValueError(f"Unsupported station: {station}")
    if any(t.stats.station != station for t in stream):
        raise ValueError("Input contains a different station")
    from tiskitpy import CleanRotator, PeriodicTransient

    period, clips = PARAMETERS[station]
    output = Stream()
    for chunk in split_stream(stream, window_seconds):
        if any(t.stats.npts / t.stats.sampling_rate < 2 * period for t in chunk):
            raise ValueError(
                "Each chunk needs at least two transient periods; adjust window size"
            )
        if rotate:
            chunk = CleanRotator(
                chunk, remove_eqs=earthquake_spans, filt_band=(0.001, 1)
            ).apply(chunk, horiz_too=True)
        transient = PeriodicTransient("hourly", period, 0.05, clips, transient_start)

        def cleaner(trace, transient=transient):
            if interactive:
                transient.calc_timing(trace, eq_spans=earthquake_spans)
            transient.calc_transient(trace, eq_spans=earthquake_spans, plots=False)
            return transient.remove_transient(trace, match=False, plots=False)

        output += clean_chunk(chunk, cleaner)
    output.merge(method=0)
    output.detrend("demean")
    output.detrend("linear")
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--station", choices=sorted(PARAMETERS), required=True)
    parser.add_argument(
        "--transient-start", required=True, help="Calibrated UTC transient onset"
    )
    parser.add_argument("--window-days", type=float, default=5)
    parser.add_argument("--rotate", action="store_true")
    parser.add_argument("--interactive", action="store_true")
    parser.add_argument("--earthquake-spans", action="store_true")
    args = parser.parse_args(argv)
    if args.output.exists():
        parser.error("output already exists; choose a new path")
    cleaned = process(
        read(str(args.input)),
        args.station,
        UTCDateTime(args.transient_start),
        window_seconds=args.window_days * 86400,
        rotate=args.rotate,
        interactive=args.interactive,
        earthquake_spans=args.earthquake_spans,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for trace in cleaned:
        trace.data = np.asarray(trace.data, dtype=np.float64)
    cleaned.write(str(args.output), format="MSEED")


if __name__ == "__main__":
    main()
