import unittest

try:
    import importlib.resources as importlib_resources
except ImportError:
    import importlib_resources

import pymuscle


class TestAligner(unittest.TestCase):

    def test_halorhodopsin(self):

        with importlib_resources.path("pymuscle.tests.data", "swissprot-halorhodopsin.faa") as path:
            sequences = pymuscle.MultiSequence.from_file(path)
            self.assertEqual(len(sequences), 11)

        aligner = pymuscle.Aligner()
        msa = aligner.align(sequences)
        self.assertEqual(len(msa.sequences), 11)

        with importlib_resources.path("pymuscle.tests.data", "swissprot-halorhodopsin.muscle.afa") as path:
            expected = pymuscle.MultiSequence.from_file(path)
            self.assertEqual(len(expected), 11)

        actual_sequences = { seq.name:seq for seq in msa.sequences }
        expected_sequences = {seq.name:seq for seq in expected }

        self.assertEqual(set(actual_sequences), set(expected_sequences))
        for name in sorted(expected_sequences):
            actual_seq = actual_sequences[name]
            expected_seq = expected_sequences[name]
            self.assertEqual(actual_seq.name.decode(), expected_seq.name.decode())
            self.assertEqual(actual_seq.sequence.decode(), expected_seq.sequence.decode())
