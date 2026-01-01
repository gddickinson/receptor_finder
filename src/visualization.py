"""
Visualization Module
====================

Plotting and visualization tools for receptor finder analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Optional
import pandas as pd
from pathlib import Path


class CandidatePlotter:
    """Plotting tools for candidate analysis."""
    
    @staticmethod
    def plot_score_distribution(scores: pd.DataFrame,
                               output_path: Optional[str] = None,
                               title: str = "Candidate Score Distribution"):
        """
        Plot distribution of candidate scores.
        
        Parameters
        ----------
        scores : pd.DataFrame
            DataFrame with candidate scores
        output_path : Optional[str]
            Save path
        title : str
            Plot title
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Total score distribution
        axes[0, 0].hist(scores['total_score'], bins=20, color='steelblue', 
                       edgecolor='black', alpha=0.7)
        axes[0, 0].axvline(60, color='red', linestyle='--', linewidth=2, 
                          label='High Priority Threshold')
        axes[0, 0].axvline(40, color='orange', linestyle='--', linewidth=2,
                          label='Medium Priority Threshold')
        axes[0, 0].set_xlabel('Total Score', fontsize=12, fontweight='bold')
        axes[0, 0].set_ylabel('Number of Candidates', fontsize=12, fontweight='bold')
        axes[0, 0].set_title('Total Score Distribution', fontweight='bold')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Priority class pie chart
        priority_counts = scores['priority_class'].value_counts()
        colors = {'HIGH': 'green', 'MEDIUM': 'orange', 'LOW': 'gray'}
        pie_colors = [colors.get(p, 'gray') for p in priority_counts.index]
        axes[0, 1].pie(priority_counts.values, labels=priority_counts.index,
                      autopct='%1.1f%%', colors=pie_colors, startangle=90)
        axes[0, 1].set_title('Priority Classification', fontweight='bold')
        
        # Score component comparison
        score_cols = ['pocket_score', 'docking_score', 'structure_quality',
                     'phylo_pattern_score', 'coexpression_score', 'localization_score']
        score_means = scores[score_cols].mean()
        
        axes[1, 0].barh(range(len(score_means)), score_means.values, color='teal', alpha=0.7)
        axes[1, 0].set_yticks(range(len(score_means)))
        axes[1, 0].set_yticklabels([c.replace('_score', '').replace('_', ' ').title() 
                                    for c in score_means.index])
        axes[1, 0].set_xlabel('Average Score', fontsize=12, fontweight='bold')
        axes[1, 0].set_title('Score Components', fontweight='bold')
        axes[1, 0].grid(True, alpha=0.3, axis='x')
        
        # Top candidates
        top_candidates = scores.nlargest(10, 'total_score')
        axes[1, 1].barh(range(len(top_candidates)), top_candidates['total_score'].values,
                       color='darkgreen', alpha=0.7)
        axes[1, 1].set_yticks(range(len(top_candidates)))
        axes[1, 1].set_yticklabels(top_candidates['protein_id'].values, fontsize=8)
        axes[1, 1].set_xlabel('Total Score', fontsize=12, fontweight='bold')
        axes[1, 1].set_title('Top 10 Candidates', fontweight='bold')
        axes[1, 1].grid(True, alpha=0.3, axis='x')
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
    
    @staticmethod
    def plot_score_heatmap(scores: pd.DataFrame,
                          output_path: Optional[str] = None,
                          title: str = "Candidate Score Heatmap"):
        """
        Plot heatmap of score components for top candidates.
        
        Parameters
        ----------
        scores : pd.DataFrame
            DataFrame with candidate scores
        output_path : Optional[str]
            Save path
        title : str
            Plot title
        """
        # Select top 20 candidates
        top_candidates = scores.nlargest(20, 'total_score')
        
        # Select score columns
        score_cols = ['pocket_score', 'docking_score', 'structure_quality',
                     'phylo_pattern_score', 'expansion_score',
                     'coexpression_score', 'localization_score', 'tm_domain_score']
        
        # Create matrix
        matrix = top_candidates[score_cols].values
        
        # Normalize by max possible score for each component
        max_scores = np.array([13, 10, 7, 10, 10, 10, 10, 5])
        matrix_norm = matrix / max_scores
        
        # Plot
        fig, ax = plt.subplots(figsize=(10, 12))
        
        sns.heatmap(matrix_norm, annot=False, cmap='RdYlGn', vmin=0, vmax=1,
                   xticklabels=[c.replace('_score', '').replace('_', ' ').title() 
                               for c in score_cols],
                   yticklabels=top_candidates['protein_id'].values,
                   cbar_kws={'label': 'Normalized Score (0-1)'},
                   ax=ax)
        
        ax.set_xlabel('Score Component', fontsize=12, fontweight='bold')
        ax.set_ylabel('Protein ID', fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold')
        
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
    
    @staticmethod
    def plot_docking_comparison(docking_results: pd.DataFrame,
                               output_path: Optional[str] = None,
                               title: str = "Ligand Docking Comparison"):
        """
        Plot comparison of docking results for different ligands.
        
        Parameters
        ----------
        docking_results : pd.DataFrame
            DataFrame with docking results
        output_path : Optional[str]
            Save path
        title : str
            Plot title
        """
        # Group by protein and ligand
        pivot_data = docking_results.pivot_table(
            values='binding_affinity',
            index='protein_id',
            columns='ligand_id',
            aggfunc='min'  # Best affinity
        )
        
        # Select top 15 proteins by best InsP3 affinity
        if 'InsP3' in pivot_data.columns:
            top_proteins = pivot_data.nsmallest(15, 'InsP3')
        else:
            top_proteins = pivot_data.head(15)
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        x = np.arange(len(top_proteins))
        width = 0.2
        
        ligands = top_proteins.columns
        colors = ['green', 'blue', 'orange', 'red', 'purple', 'brown']
        
        for i, ligand in enumerate(ligands):
            offset = (i - len(ligands)/2) * width
            ax.bar(x + offset, top_proteins[ligand].values, width,
                  label=ligand, color=colors[i % len(colors)], alpha=0.7)
        
        ax.axhline(-7, color='red', linestyle='--', linewidth=2,
                  label='Good Binding Threshold')
        
        ax.set_xlabel('Protein ID', fontsize=12, fontweight='bold')
        ax.set_ylabel('Binding Affinity (kcal/mol)', fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(top_proteins.index, rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
    
    @staticmethod
    def plot_structure_quality(structures: pd.DataFrame,
                              output_path: Optional[str] = None,
                              title: str = "Structure Quality Assessment"):
        """
        Plot structure quality metrics.
        
        Parameters
        ----------
        structures : pd.DataFrame
            DataFrame with structure quality data
        output_path : Optional[str]
            Save path
        title : str
            Plot title
        """
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # pLDDT distribution
        axes[0].hist(structures['mean_plddt'], bins=20, color='purple',
                    edgecolor='black', alpha=0.7)
        axes[0].axvline(70, color='red', linestyle='--', linewidth=2,
                       label='Quality Threshold')
        axes[0].axvline(90, color='green', linestyle='--', linewidth=2,
                       label='High Confidence')
        axes[0].set_xlabel('Mean pLDDT', fontsize=12, fontweight='bold')
        axes[0].set_ylabel('Number of Structures', fontsize=12, fontweight='bold')
        axes[0].set_title('AlphaFold2 Confidence Scores', fontweight='bold')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # Quality vs TM quality scatter
        if 'tm_region_quality' in structures.columns:
            valid = structures['tm_region_quality'].notna()
            axes[1].scatter(structures.loc[valid, 'mean_plddt'],
                          structures.loc[valid, 'tm_region_quality'],
                          c=structures.loc[valid, 'is_high_quality'].map({True: 'green', False: 'red'}),
                          alpha=0.6, s=50, edgecolors='black')
            axes[1].plot([0, 100], [0, 100], 'k--', alpha=0.3)
            axes[1].set_xlabel('Mean pLDDT', fontsize=12, fontweight='bold')
            axes[1].set_ylabel('TM Region pLDDT', fontsize=12, fontweight='bold')
            axes[1].set_title('Overall vs TM Region Quality', fontweight='bold')
            axes[1].grid(True, alpha=0.3)
            axes[1].set_xlim(0, 100)
            axes[1].set_ylim(0, 100)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
    
    @staticmethod
    def plot_pocket_characteristics(pockets: pd.DataFrame,
                                   output_path: Optional[str] = None,
                                   title: str = "Binding Pocket Characteristics"):
        """
        Plot binding pocket characteristics.
        
        Parameters
        ----------
        pockets : pd.DataFrame
            DataFrame with pocket data
        output_path : Optional[str]
            Save path
        title : str
            Plot title
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Volume vs charge
        axes[0, 0].scatter(pockets['volume'], pockets['charge'],
                          c=pockets['insp3_score'], cmap='viridis',
                          s=100, alpha=0.6, edgecolors='black')
        axes[0, 0].set_xlabel('Pocket Volume (ų)', fontsize=12, fontweight='bold')
        axes[0, 0].set_ylabel('Net Charge', fontsize=12, fontweight='bold')
        axes[0, 0].set_title('Volume vs Charge', fontweight='bold')
        axes[0, 0].grid(True, alpha=0.3)
        cbar = plt.colorbar(axes[0, 0].collections[0], ax=axes[0, 0])
        cbar.set_label('InsP₃ Score', fontsize=10)
        
        # InsP3 vs cADPR scores
        axes[0, 1].scatter(pockets['insp3_score'], pockets['cadpr_score'],
                          alpha=0.6, s=100, edgecolors='black', color='teal')
        axes[0, 1].plot([0, 13], [0, 13], 'k--', alpha=0.3)
        axes[0, 1].set_xlabel('InsP₃ Score', fontsize=12, fontweight='bold')
        axes[0, 1].set_ylabel('cADPR Score', fontsize=12, fontweight='bold')
        axes[0, 1].set_title('Ligand Binding Scores', fontweight='bold')
        axes[0, 1].grid(True, alpha=0.3)
        
        # Residue composition
        residue_types = ['basic_count', 'aromatic_count', 'acidic_count', 'polar_count']
        means = pockets[residue_types].mean()
        
        axes[1, 0].bar(range(len(means)), means.values, color='coral', alpha=0.7,
                      edgecolor='black')
        axes[1, 0].set_xticks(range(len(means)))
        axes[1, 0].set_xticklabels([r.replace('_count', '').title() for r in residue_types])
        axes[1, 0].set_ylabel('Average Count', fontsize=12, fontweight='bold')
        axes[1, 0].set_title('Pocket Residue Composition', fontweight='bold')
        axes[1, 0].grid(True, alpha=0.3, axis='y')
        
        # Score distribution
        axes[1, 1].hist(pockets['insp3_score'], bins=15, alpha=0.6, 
                       label='InsP₃', color='green', edgecolor='black')
        axes[1, 1].hist(pockets['cadpr_score'], bins=15, alpha=0.6,
                       label='cADPR', color='blue', edgecolor='black')
        axes[1, 1].set_xlabel('Binding Score', fontsize=12, fontweight='bold')
        axes[1, 1].set_ylabel('Number of Pockets', fontsize=12, fontweight='bold')
        axes[1, 1].set_title('Binding Score Distribution', fontweight='bold')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            plt.close()
        else:
            plt.show()


if __name__ == "__main__":
    # Test code
    print("Testing Visualization Module\n")
    
    # Create mock data
    print("Generating mock data...")
    
    # Mock candidate scores
    np.random.seed(42)
    n_candidates = 50
    
    scores_data = {
        'protein_id': [f'AT{i}G{12345+i}' for i in range(n_candidates)],
        'total_score': np.random.uniform(20, 90, n_candidates),
        'pocket_score': np.random.uniform(5, 13, n_candidates),
        'docking_score': np.random.uniform(3, 10, n_candidates),
        'structure_quality': np.random.uniform(2, 7, n_candidates),
        'phylo_pattern_score': np.random.uniform(2, 10, n_candidates),
        'expansion_score': np.random.uniform(2, 10, n_candidates),
        'coexpression_score': np.random.uniform(2, 10, n_candidates),
        'localization_score': np.random.uniform(2, 10, n_candidates),
        'tm_domain_score': np.random.uniform(1, 5, n_candidates),
    }
    
    scores_df = pd.DataFrame(scores_data)
    scores_df['priority_class'] = scores_df['total_score'].apply(
        lambda x: 'HIGH' if x >= 60 else ('MEDIUM' if x >= 40 else 'LOW')
    )
    
    # Test 1: Score distribution
    print("\nTest 1: Score Distribution Plot")
    CandidatePlotter.plot_score_distribution(
        scores_df,
        output_path='test_outputs/score_distribution.png'
    )
    print("  Saved to test_outputs/score_distribution.png")
    
    # Test 2: Score heatmap
    print("\nTest 2: Score Heatmap")
    CandidatePlotter.plot_score_heatmap(
        scores_df,
        output_path='test_outputs/score_heatmap.png'
    )
    print("  Saved to test_outputs/score_heatmap.png")
    
    print("\nVisualization tests complete!")
