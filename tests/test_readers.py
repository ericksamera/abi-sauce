from abi_sauce.io import readers

FASTA = b">x\nACGT\n"

def test_fasta_single():
    out = readers.read_fasta("x.fa", FASTA, len(FASTA), "md5")
    assert len(out) == 1
    assert out[0].sequence == "ACGT"
