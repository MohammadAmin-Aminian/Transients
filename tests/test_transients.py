import numpy as np
import pytest
from obspy import Stream, Trace
from Transients_Removal_24_Aug_23 import clean_chunk, split_stream


def stream():
    return Stream(
        [
            Trace(
                np.arange(11, dtype=float),
                header={"channel": channel, "sampling_rate": 1},
            )
            for channel in ["BHZ", "BH1"]
        ]
    )


def test_partition_preserves_every_sample_once():
    chunks = list(split_stream(stream(), 4))
    assert [c[0].stats.npts for c in chunks] == [4, 4, 3]
    assert [c[0].stats.starttime.timestamp for c in chunks] == [0, 4, 8]
    merged = sum(chunks, Stream()).merge()
    np.testing.assert_array_equal(merged.select(channel="BHZ")[0].data, np.arange(11))


def test_cleaner_replaces_parent_not_selection():
    original = stream()

    def cleaner(trace):
        trace.data[:] = 0
        trace.stats.channel = "CLN"
        return trace

    result = clean_chunk(original, cleaner)
    assert not result.select(channel="BHZ")[0].data.any()
    np.testing.assert_array_equal(result.select(channel="BH1")[0].data, np.arange(11))
    assert original[0].data.sum() > 0


def test_misalignment_rejected():
    data = stream()
    data[1].stats.starttime += 1
    with pytest.raises(ValueError):
        list(split_stream(data, 4))


def test_process_appends_each_cleaned_chunk_once(monkeypatch):
    import tiskitpy
    from Transients_Removal_24_Aug_23 import process

    class Transient:
        def __init__(self, *args):
            pass

        def calc_transient(self, trace, **kwargs):
            pass

        def remove_transient(self, trace, **kwargs):
            trace.data *= 0
            return trace

    monkeypatch.setattr(tiskitpy, "PeriodicTransient", Transient)
    data = Stream(
        [
            Trace(
                np.ones(16000),
                header={"station": "RR38", "channel": "BHZ", "sampling_rate": 1},
            )
        ]
    )
    output = process(data, "RR38", data[0].stats.starttime, window_seconds=8000)
    assert len(output) == 1
    assert output[0].stats.npts == 16000
    assert not output[0].data.any()


def test_real_tiskitpy_suppresses_known_periodic_pulses():
    from obspy import UTCDateTime
    from Transients_Removal_24_Aug_23 import process

    samples = 16000
    time = np.arange(samples)
    background = np.random.default_rng(2).normal(0, 0.1, samples)
    pulses = np.zeros(samples)
    for onset in range(1000, samples, 3620):
        pulses += 10 * np.exp(-0.5 * ((time - onset) / 15) ** 2)
    original = background + pulses
    data = Stream(
        [
            Trace(
                original.copy(),
                header={"station": "RR38", "channel": "BHZ", "sampling_rate": 1},
            )
        ]
    )
    cleaned = process(data, "RR38", UTCDateTime(1000), window_seconds=samples)
    assert len(cleaned) == 1
    assert cleaned[0].id == data[0].id
    for field in ("starttime", "endtime", "sampling_rate", "npts"):
        assert cleaned[0].stats[field] == data[0].stats[field]
    assert np.isfinite(cleaned[0].data).all()
    before = np.sqrt(np.mean(pulses**2))
    after = np.sqrt(np.mean((cleaned[0].data - background) ** 2))
    assert after < before * 0.2
    np.testing.assert_array_equal(data[0].data, original)
