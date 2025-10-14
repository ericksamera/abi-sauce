from abi_sauce.services.file_manager import FileManager


def test_checksum_dedup():
    fm = FileManager()
    fasta = b">x\nACGT\n"
    out1 = fm.add_bytes(name="x.fa", raw=fasta)
    out2 = fm.add_bytes(name="copy.fa", raw=fasta)
    assert len(out1) == 1
    assert out2 == []  # duplicate bytes skipped
