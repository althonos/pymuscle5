import unittest

try:
    import importlib.resources as importlib_resources
except ImportError:
    import importlib_resources

from .. import Aligner, MultiSequence
from .data import __name__ as data_package


class TestAligner(unittest.TestCase):

    def test_luxc(self):

        with importlib_resources.path(data_package, "LuxC.faa") as path:
            sequences = MultiSequence.from_file(str(path))
            self.assertEqual(len(sequences), 13)

        aligner = Aligner()
        msa = aligner.align(sequences)
        self.assertEqual(len(msa.sequences), 13)

        with importlib_resources.path(data_package, "LuxC.muscle.afa") as path:
            expected = MultiSequence.from_file(str(path))
            self.assertEqual(len(expected), 13)

        actual_sequences = { seq.name:seq for seq in msa.sequences }
        expected_sequences = {seq.name:seq for seq in expected }

        self.assertEqual(set(actual_sequences), set(expected_sequences))
        for name in sorted(expected_sequences):
            actual_seq = actual_sequences[name]
            expected_seq = expected_sequences[name]
            self.assertEqual(actual_seq.name.decode(), expected_seq.name.decode())
            self.assertEqual(actual_seq.sequence.decode(), expected_seq.sequence.decode())

    def test_halorhodopsin(self):

        with importlib_resources.path(data_package, "swissprot-halorhodopsin.faa") as path:
            sequences = MultiSequence.from_file(str(path))
            self.assertEqual(len(sequences), 11)

        aligner = Aligner()
        msa = aligner.align(sequences)
        self.assertEqual(len(msa.sequences), 11)

        with importlib_resources.path(data_package, "swissprot-halorhodopsin.muscle.afa") as path:
            expected = MultiSequence.from_file(str(path))
            self.assertEqual(len(expected), 11)

        actual_sequences = { seq.name:seq for seq in msa.sequences }
        expected_sequences = {seq.name:seq for seq in expected }

        self.assertEqual(set(actual_sequences), set(expected_sequences))
        for name in sorted(expected_sequences):
            actual_seq = actual_sequences[name]
            expected_seq = expected_sequences[name]
            self.assertEqual(actual_seq.name.decode(), expected_seq.name.decode())
            self.assertEqual(actual_seq.sequence.decode(), expected_seq.sequence.decode())
