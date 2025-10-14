from abi_sauce.models import Sample
from abi_sauce.services.file_manager import FileManager
from abi_sauce.services.sample_manager import SampleManager


def test_sequence_override_semantics():
    fm = FileManager()
    sm = SampleManager(fm)
    s = Sample(id="1", name="x")
    sm._samples["1"] = s

    sm.set_sequence_override("1", None)
    assert sm.get("1").sequence_override is None
    sm.set_sequence_override("1", "")
    assert sm.get("1").sequence_override is None
    sm.set_sequence_override("1", "AC\nGT\r")
    assert sm.get("1").sequence_override == "ACGT"


def test_tags_editing():
    fm = FileManager()
    sm = SampleManager(fm)
    s = Sample(id="1", name="x")
    sm._samples["1"] = s

    sm.set_tags("1", ["  foo ", "", "bar"])
    assert sm.get("1").tags == ["foo", "bar"]
