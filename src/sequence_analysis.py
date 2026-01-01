"""
Sequence Analysis Module
=========================

Handles protein sequence analysis including:
- Sequence retrieval and filtering
- Transmembrane domain prediction
- BLAST-like homology searches
- HMM profile matching
- Candidate protein identification
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import re
from pathlib import Path


@dataclass
class ProteinSequence:
    """Data class for protein sequences."""
    
    sequence_id: str
    species: str
    sequence: str
    description: Optional[str] = None
    
    # Predicted features
    n_transmembrane: Optional[int] = None
    localization: Optional[str] = None
    has_signal_peptide: bool = False
    
    # Annotations
    domains: Optional[List[str]] = None
    gene_family: Optional[str] = None
    
    def __post_init__(self):
        """Calculate basic properties."""
        self.length = len(self.sequence)
        
    @property
    def molecular_weight(self) -> float:
        """Calculate approximate molecular weight (kDa)."""
        # Average residue weight ~110 Da
        return self.length * 0.11
    
    def count_residue(self, residue: str) -> int:
        """Count specific residue."""
        return self.sequence.count(residue.upper())
    
    def has_motif(self, motif: str) -> bool:
        """Check if sequence contains motif."""
        return motif.upper() in self.sequence.upper()


class TransmembraneDomainPredictor:
    """Simple transmembrane domain prediction based on hydrophobicity."""
    
    # Kyte-Doolittle hydrophobicity scale
    HYDROPHOBICITY = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }
    
    @staticmethod
    def predict_tm_domains(sequence: str, window: int = 19, 
                          threshold: float = 1.6) -> List[Tuple[int, int]]:
        """
        Predict transmembrane domains using sliding window.
        
        Parameters
        ----------
        sequence : str
            Protein sequence
        window : int
            Window size for hydrophobicity calculation
        threshold : float
            Hydrophobicity threshold for TM domain
            
        Returns
        -------
        tm_domains : List[Tuple[int, int]]
            List of (start, end) positions for TM domains
        """
        seq = sequence.upper()
        n = len(seq)
        
        if n < window:
            return []
        
        # Calculate hydrophobicity profile
        hydro_profile = []
        for i in range(n - window + 1):
            window_seq = seq[i:i+window]
            avg_hydro = np.mean([
                TransmembraneDomainPredictor.HYDROPHOBICITY.get(aa, 0) 
                for aa in window_seq
            ])
            hydro_profile.append(avg_hydro)
        
        # Find regions above threshold
        tm_domains = []
        in_tm = False
        start = 0
        
        for i, hydro in enumerate(hydro_profile):
            if hydro >= threshold and not in_tm:
                start = i
                in_tm = True
            elif hydro < threshold and in_tm:
                # Require minimum length of 15 residues
                if i - start >= 15:
                    tm_domains.append((start, i + window))
                in_tm = False
        
        # Check if still in TM at end
        if in_tm and len(hydro_profile) - start >= 15:
            tm_domains.append((start, len(hydro_profile) + window))
        
        # Merge nearby TM domains (within 5 residues)
        if len(tm_domains) > 1:
            merged = [tm_domains[0]]
            for current in tm_domains[1:]:
                last = merged[-1]
                if current[0] - last[1] <= 5:
                    merged[-1] = (last[0], current[1])
                else:
                    merged.append(current)
            tm_domains = merged
        
        return tm_domains


class SequenceFilter:
    """Filter protein sequences based on criteria."""
    
    @staticmethod
    def filter_by_length(sequences: List[ProteinSequence],
                        min_length: int = 200,
                        max_length: int = 2000) -> List[ProteinSequence]:
        """Filter sequences by length."""
        return [seq for seq in sequences 
                if min_length <= seq.length <= max_length]
    
    @staticmethod
    def filter_by_tm_domains(sequences: List[ProteinSequence],
                            min_tm: int = 2,
                            max_tm: int = 10) -> List[ProteinSequence]:
        """Filter sequences by number of TM domains."""
        filtered = []
        predictor = TransmembraneDomainPredictor()
        
        for seq in sequences:
            tm_domains = predictor.predict_tm_domains(seq.sequence)
            n_tm = len(tm_domains)
            seq.n_transmembrane = n_tm
            
            if min_tm <= n_tm <= max_tm:
                filtered.append(seq)
        
        return filtered
    
    @staticmethod
    def filter_by_motifs(sequences: List[ProteinSequence],
                        required_motifs: Optional[List[str]] = None,
                        excluded_motifs: Optional[List[str]] = None) -> List[ProteinSequence]:
        """Filter sequences by presence/absence of motifs."""
        filtered = []
        
        for seq in sequences:
            # Check required motifs
            if required_motifs:
                has_all = all(seq.has_motif(motif) for motif in required_motifs)
                if not has_all:
                    continue
            
            # Check excluded motifs
            if excluded_motifs:
                has_any = any(seq.has_motif(motif) for motif in excluded_motifs)
                if has_any:
                    continue
            
            filtered.append(seq)
        
        return filtered


class MotifScanner:
    """Scan for specific sequence motifs."""
    
    # Common motifs for Ca2+ channels and InsP3/cADPR binding
    MOTIFS = {
        # Nucleotide binding motifs
        'P_loop': r'G[KR].{1,5}G[KR][ST]',  # Walker A
        'Walker_B': r'[RK].{2,4}[LIVM].{2,4}D',
        
        # Ca2+ binding
        'EF_hand': r'D.{2}[DNS]...[DNS].{2}[EQ]',
        'Acidic_cluster': r'[DE]{3,5}',  # Multiple acidic residues
        
        # Phosphate binding
        'Arg_cluster': r'[RK]{2,4}',  # Basic cluster for phosphates
        
        # Pore motifs
        'Selectivity_filter': r'[TSD][TSD][DE]',
    }
    
    @staticmethod
    def scan_sequence(sequence: str, motif_pattern: str) -> List[int]:
        """
        Scan sequence for motif pattern.
        
        Returns list of start positions.
        """
        positions = []
        for match in re.finditer(motif_pattern, sequence, re.IGNORECASE):
            positions.append(match.start())
        return positions
    
    @staticmethod
    def get_motif_counts(sequence: str) -> Dict[str, int]:
        """Count occurrences of all known motifs."""
        counts = {}
        for motif_name, pattern in MotifScanner.MOTIFS.items():
            matches = MotifScanner.scan_sequence(sequence, pattern)
            counts[motif_name] = len(matches)
        return counts


class SequenceDatabase:
    """Manage protein sequence database."""
    
    def __init__(self):
        self.sequences: List[ProteinSequence] = []
    
    def add_sequence(self, sequence: ProteinSequence):
        """Add sequence to database."""
        self.sequences.append(sequence)
    
    def load_from_fasta(self, filepath: str, species: str = "Unknown"):
        """
        Load sequences from FASTA file.
        
        Parameters
        ----------
        filepath : str
            Path to FASTA file
        species : str
            Species name for all sequences
        """
        current_id = None
        current_desc = None
        current_seq = []
        
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence
                    if current_id:
                        seq = ProteinSequence(
                            sequence_id=current_id,
                            species=species,
                            sequence=''.join(current_seq),
                            description=current_desc
                        )
                        self.add_sequence(seq)
                    
                    # Parse header
                    header = line[1:]
                    parts = header.split(None, 1)
                    current_id = parts[0]
                    current_desc = parts[1] if len(parts) > 1 else None
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Save last sequence
            if current_id:
                seq = ProteinSequence(
                    sequence_id=current_id,
                    species=species,
                    sequence=''.join(current_seq),
                    description=current_desc
                )
                self.add_sequence(seq)
    
    def save_to_fasta(self, filepath: str):
        """Save sequences to FASTA file."""
        with open(filepath, 'w') as f:
            for seq in self.sequences:
                header = f">{seq.sequence_id}"
                if seq.description:
                    header += f" {seq.description}"
                f.write(header + "\n")
                
                # Write sequence in 60-character lines
                for i in range(0, len(seq.sequence), 60):
                    f.write(seq.sequence[i:i+60] + "\n")
    
    def to_dataframe(self) -> pd.DataFrame:
        """Convert to pandas DataFrame."""
        data = []
        for seq in self.sequences:
            data.append({
                'sequence_id': seq.sequence_id,
                'species': seq.species,
                'length': seq.length,
                'n_transmembrane': seq.n_transmembrane,
                'localization': seq.localization,
                'gene_family': seq.gene_family,
                'description': seq.description
            })
        return pd.DataFrame(data)
    
    def filter_database(self, 
                       min_length: int = 200,
                       max_length: int = 2000,
                       min_tm: int = 2,
                       max_tm: int = 10) -> 'SequenceDatabase':
        """
        Filter database with multiple criteria.
        
        Returns new SequenceDatabase with filtered sequences.
        """
        filtered_db = SequenceDatabase()
        
        # Apply filters sequentially
        seqs = self.sequences
        seqs = SequenceFilter.filter_by_length(seqs, min_length, max_length)
        seqs = SequenceFilter.filter_by_tm_domains(seqs, min_tm, max_tm)
        
        filtered_db.sequences = seqs
        return filtered_db
    
    def __len__(self):
        return len(self.sequences)
    
    def __getitem__(self, idx):
        return self.sequences[idx]


def create_example_database() -> SequenceDatabase:
    """Create example database with hypothetical candidates."""
    db = SequenceDatabase()
    
    # Example 1: Hypothetical TM protein with basic clusters
    seq1 = ProteinSequence(
        sequence_id="AT1G12345",
        species="Arabidopsis thaliana",
        sequence="MAASLKRRKRGSTTAAALLLVVVVGGGIIIIFFFLLLAAARRRKKK" +
                 "DEEEDDDGGSSHHHTTTNNNQQQLLLIIIFFFVVVAAAGGGRRRKKKK" +
                 "EEEDDDSSSTTTNNNGGGLLLIIIFFFVVVAAASSSDDDEEEKKKRRR" +
                 "GGGLLLIIIFFFVVVAAATTTSSSDDDEEEGGGRRRKKKHHHTTTAAA" +
                 "LLLIIIFFVVVGGGAAASSSTTTEEEDDDKKRRRGGGHHHTTTNNNAAA" +
                 "LLLIIFFFVVGGGAAASSSDEEEKKKRRRGGGHHHTTTNNQQQ",
        description="Hypothetical membrane protein with basic clusters",
        gene_family="Unknown"
    )
    
    # Example 2: Protein with EF-hand motifs
    seq2 = ProteinSequence(
        sequence_id="AT2G23456",
        species="Arabidopsis thaliana",
        sequence="MAAADKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTI" +
                 "DFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTN" +
                 "LGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAKGGGKKKRRRAAA" +
                 "LLLIIIFFFVVVGGGSSSDDDEEELLLIIFFFVVAAAGGGKKRRRHHH",
        description="Protein with Ca2+ binding motifs",
        gene_family="EF-hand"
    )
    
    # Example 3: Multi-TM protein with phosphate binding potential
    seq3 = ProteinSequence(
        sequence_id="AT3G34567",
        species="Arabidopsis thaliana",
        sequence="MAAAGKSGKSTRRRKKKVVVLLLIIIFFAAALLLIIIFFVVVAAAGGG" +
                 "SSSDDDEEELLLIIIFFVVVAAAGGGKKRRRSSSDDDEEELLLIIIFFF" +
                 "VVVAAAGGGKKRRRHHHTTTNNNQQQEEEDDDLLLIIIFFVVVAAAGGG" +
                 "KKRRRSSSDDDEEEGGGLLLIIIFFVVVAAATTTSSSGGGKKRRREEEDD",
        description="Multi-pass TM protein",
        gene_family="Unknown_TM"
    )
    
    for seq in [seq1, seq2, seq3]:
        # Predict TM domains
        predictor = TransmembraneDomainPredictor()
        tm_domains = predictor.predict_tm_domains(seq.sequence)
        seq.n_transmembrane = len(tm_domains)
        
        db.add_sequence(seq)
    
    return db


if __name__ == "__main__":
    # Test code
    print("Testing Sequence Analysis Module\n")
    
    # Test 1: TM domain prediction
    print("Test 1: TM Domain Prediction")
    test_seq = "MAAALLLIIIFFFVVVGGGAAASSSLLLIIIFFVVVGGGAAASSS"
    predictor = TransmembraneDomainPredictor()
    tm_domains = predictor.predict_tm_domains(test_seq)
    print(f"  Sequence length: {len(test_seq)}")
    print(f"  TM domains found: {len(tm_domains)}")
    for i, (start, end) in enumerate(tm_domains):
        print(f"    TM{i+1}: {start}-{end}")
    
    # Test 2: Motif scanning
    print("\nTest 2: Motif Scanning")
    test_seq2 = "MAAAGKSGKSTRRRKKKDEEEDDDGGGKKK"
    motif_counts = MotifScanner.get_motif_counts(test_seq2)
    print(f"  Sequence: {test_seq2}")
    print(f"  Motifs found:")
    for motif, count in motif_counts.items():
        if count > 0:
            print(f"    {motif}: {count}")
    
    # Test 3: Example database
    print("\nTest 3: Example Database")
    db = create_example_database()
    print(f"  Created database with {len(db)} sequences")
    
    df = db.to_dataframe()
    print("\n  Summary:")
    print(df[['sequence_id', 'length', 'n_transmembrane', 'gene_family']])
    
    # Test 4: Filtering
    print("\nTest 4: Filtering")
    filtered_db = db.filter_database(min_tm=2, max_tm=6)
    print(f"  Sequences with 2-6 TM domains: {len(filtered_db)}")
    
    # Save example
    output_dir = Path(__file__).parent.parent / "data"
    output_dir.mkdir(exist_ok=True)
    db.save_to_fasta(output_dir / "example_sequences.fasta")
    print(f"\n  Saved to {output_dir / 'example_sequences.fasta'}")
