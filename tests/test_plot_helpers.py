# tests/test_plot_helpers.py
from abi_sauce.services.trace_processing import (
    compute_base_windows,
    resample_channel_to_bases,
)


def test_compute_base_windows_simple():
    ploc = [10, 20, 30]
    total_samples = 40
    windows = compute_base_windows(ploc, total_samples)
    assert len(windows) == 3
    for s, e in windows:
        assert 0 <= s <= e < total_samples
    assert windows[0][0] == 0
    assert windows[-1][1] == total_samples - 1


def test_resample_channel_to_bases_max_mean():
    channel = list(range(40))  # 0..39
    ploc = [10, 20, 30]
    windows = compute_base_windows(ploc, total_samples=40)
    res_max = resample_channel_to_bases(channel, windows, method="max")
    res_mean = resample_channel_to_bases(channel, windows, method="mean")
    # Using midpoint partitioning, maxima should be near mid of each window
    assert 14 <= res_max[0] <= 16
    assert 24 <= res_max[1] <= 26
    assert res_max[2] == 39
    for m, M in zip(res_mean, res_max, strict=False):
        assert m <= M
