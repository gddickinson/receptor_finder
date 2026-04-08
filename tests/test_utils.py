"""
Unit tests for the utils module.
"""

import sys
import os
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from utils import read_fasta, write_fasta, validate_fasta, validate_pdb


class TestFastaIO(unittest.TestCase):
    def test_round_trip(self):
        seqs = {"SEQ1": "MAAALL", "SEQ2": "GGGKKK"}
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as f:
            tmp = f.name
        try:
            write_fasta(seqs, tmp)
            loaded = read_fasta(tmp)
            self.assertEqual(loaded, seqs)
        finally:
            os.unlink(tmp)


class TestValidation(unittest.TestCase):
    def test_validate_fasta_good(self):
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w") as f:
            f.write(">SEQ1\nMAAA\n")
            tmp = f.name
        try:
            ok, msg = validate_fasta(tmp)
            self.assertTrue(ok)
        finally:
            os.unlink(tmp)

    def test_validate_fasta_bad(self):
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w") as f:
            f.write("NOT_FASTA\n")
            tmp = f.name
        try:
            ok, _ = validate_fasta(tmp)
            self.assertFalse(ok)
        finally:
            os.unlink(tmp)

    def test_validate_fasta_missing(self):
        ok, _ = validate_fasta("/nonexistent/path.fasta")
        self.assertFalse(ok)

    def test_validate_pdb_missing(self):
        ok, _ = validate_pdb("/nonexistent/path.pdb")
        self.assertFalse(ok)


if __name__ == "__main__":
    unittest.main()
