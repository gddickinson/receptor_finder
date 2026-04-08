"""
Utilities Module
================

Shared helper functions for FASTA I/O, PDB parsing, and file validation.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re


# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------

def read_fasta(filepath: str) -> Dict[str, str]:
    """
    Read a FASTA file and return a dict of {header: sequence}.

    Parameters
    ----------
    filepath : str
        Path to the FASTA file.

    Returns
    -------
    sequences : Dict[str, str]
        Mapping of sequence ID to sequence string.
    """
    sequences: Dict[str, str] = {}
    current_id: Optional[str] = None
    current_seq: List[str] = []

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            sequences[current_id] = "".join(current_seq)

    return sequences


def write_fasta(sequences: Dict[str, str], filepath: str, line_width: int = 60):
    """
    Write sequences to a FASTA file.

    Parameters
    ----------
    sequences : Dict[str, str]
        Mapping of header to sequence.
    filepath : str
        Output file path.
    line_width : int
        Characters per line for sequence wrapping.
    """
    with open(filepath, "w") as fh:
        for header, seq in sequences.items():
            fh.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                fh.write(seq[i : i + line_width] + "\n")


# ---------------------------------------------------------------------------
# PDB helpers
# ---------------------------------------------------------------------------

def parse_pdb_residues(pdb_file: str) -> List[Dict]:
    """
    Extract CA atoms from a PDB file.

    Returns a list of dicts with keys: resname, chain, resid, x, y, z, bfactor.
    """
    residues = []
    with open(pdb_file, "r") as fh:
        for line in fh:
            if line.startswith("ATOM") and " CA " in line:
                residues.append({
                    "resname": line[17:20].strip(),
                    "chain": line[21],
                    "resid": int(line[22:26].strip()),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "bfactor": float(line[60:66].strip()),
                })
    return residues


def pdb_sequence(pdb_file: str) -> str:
    """
    Extract single-letter amino acid sequence from PDB CA atoms.
    """
    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    residues = parse_pdb_residues(pdb_file)
    return "".join(three_to_one.get(r["resname"], "X") for r in residues)


# ---------------------------------------------------------------------------
# File validation helpers
# ---------------------------------------------------------------------------

def validate_fasta(filepath: str) -> Tuple[bool, str]:
    """
    Check that a file looks like valid FASTA.

    Returns (is_valid, message).
    """
    path = Path(filepath)
    if not path.exists():
        return False, f"File not found: {filepath}"
    if path.stat().st_size == 0:
        return False, "File is empty"

    with open(filepath, "r") as fh:
        first_line = fh.readline().strip()
        if not first_line.startswith(">"):
            return False, "First line does not start with '>'"

    return True, "OK"


def validate_pdb(filepath: str) -> Tuple[bool, str]:
    """
    Check that a file looks like valid PDB.

    Returns (is_valid, message).
    """
    path = Path(filepath)
    if not path.exists():
        return False, f"File not found: {filepath}"
    if path.stat().st_size == 0:
        return False, "File is empty"

    with open(filepath, "r") as fh:
        has_atom = any(line.startswith("ATOM") for line in fh)

    if not has_atom:
        return False, "No ATOM records found"

    return True, "OK"


def validate_pdbqt(filepath: str) -> Tuple[bool, str]:
    """
    Check that a file looks like valid PDBQT.

    Returns (is_valid, message).
    """
    path = Path(filepath)
    if not path.exists():
        return False, f"File not found: {filepath}"
    if path.stat().st_size == 0:
        return False, "File is empty"

    return True, "OK"


if __name__ == "__main__":
    print("Testing Utils Module\n")

    # Test FASTA round-trip
    test_seqs = {"SEQ1": "MAAALLLIIIFFFVVV", "SEQ2": "GGGAAASSSLLLIIIFFFVVV"}
    tmp = "/tmp/test_utils.fasta"
    write_fasta(test_seqs, tmp)
    loaded = read_fasta(tmp)
    assert loaded == test_seqs, "FASTA round-trip failed"
    print("FASTA round-trip: OK")

    # Test validation
    ok, msg = validate_fasta(tmp)
    assert ok, msg
    print(f"FASTA validation: {msg}")

    print("\nTest complete!")
