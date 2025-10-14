from math import isnan

from abi_sauce.services.alignment_waveforms import build_aligned_waveforms


def test_build_aligned_waveforms_basic():
    # 3 columns, no gaps
    cols = [(0, 0), (1, 1), (2, 2)]
    # 4 channels, 100 samples; windows are non-overlapping thirds
    windows = [(0, 9), (10, 19), (20, 29)]
    ch = {b: list(range(30)) for b in "ACGT"}

    series, peaks = build_aligned_waveforms(
        cols, windows, ch, for_A=True, samples_per_base=5
    )
    assert len(peaks) == 3
    for base in "ACGT":
        xs, ys = series[base]
        # 5 points per column + 1 NaN separator, times 3 columns
        assert len(xs) == len(ys) == (5 + 1) * 3
        # there should be exactly 3 NaN row breaks
        assert sum(1 for v in xs if isinstance(v, float) and isnan(v)) == 3
