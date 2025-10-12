# tests/test_plot_helpers.py
from abi_sauce.ui import plot_helpers as ph


def test_compute_base_windows_simple():
    ploc = [10, 20, 30]
    total_samples = 40
    windows = ph._compute_base_windows(ploc, total_samples)
    # should produce three non-overlapping windows covering sample indices
    assert len(windows) == 3
    # ensure windows are tuples and within bounds
    for s, e in windows:
        assert 0 <= s <= e < total_samples
    # sanity checks for expected boundaries (based on midpoint splitting)
    assert windows[0][0] == 0
    assert windows[-1][1] == total_samples - 1


def test_resample_channel_to_bases_max_mean():
    # channel increasing 0..39
    channel = list(range(40))
    ploc = [10, 20, 30]
    windows = ph._compute_base_windows(ploc, total_samples=40)
    res_max = ph._resample_channel_to_bases(channel, windows, method="max")
    res_mean = ph._resample_channel_to_bases(channel, windows, method="mean")
    # For the windows produced with ploc [10,20,30], expect maxima at roughly 15,25,39
    assert (
        res_max[0] <= 15 and res_max[0] >= 14
    )  # allow 1-sample tolerance due to floor/ceil
    assert res_max[1] <= 25 and res_max[1] >= 24
    assert res_max[2] == 39
    # Mean should be less than or equal to max for each window
    for m, M in zip(res_mean, res_max):
        assert m <= M
