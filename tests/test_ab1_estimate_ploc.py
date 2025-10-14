from abi_sauce.services.ab1 import estimate_ploc


def test_estimate_ploc_basic():
    # Four obvious peaks across the sample index, one per base channel
    channels = {
        "A": [0, 0, 10, 0, 0, 0, 0],
        "C": [0, 0, 0, 0, 12, 0, 0],
        "G": [0, 0, 0, 0, 0, 9, 0],
        "T": [0, 5, 0, 0, 0, 0, 0],
    }
    ploc = estimate_ploc(channels, seq_len=4)
    assert ploc is not None
    assert len(ploc) == 4
    assert ploc == sorted(ploc)
