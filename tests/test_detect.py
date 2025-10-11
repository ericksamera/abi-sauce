from abi_sauce.io.detect import identify

def test_sniff_ab1_magic():
    kind, ext = identify("mystery.bin", b"ABIFxxxx")
    assert kind == "ab1"