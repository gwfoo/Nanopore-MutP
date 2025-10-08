#!/usr/bin/env python3
"""
Advanced Variant Analysis and Visualization for Plasmid Sequences
Supports Oxford Nanopore long-read data analysis
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam
import vcf
from collections import defaultdict, Counter
import json
import sys
import os
from pathlib import Path

class PlasmidVariantAnalyzer:
    def __init__(self, reference_file, vcf_file=None, bam_file=None, output_dir="variant_analysis"):
        self.reference_file = reference_file
        self.vcf_file = vcf_file
        self.bam_file = bam_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Load reference sequences
        self.references = SeqIO.to_dict(SeqIO.parse(reference_file, 'fasta'))
        
        # Initialize data containers
        self.variants = []
        self.structural_variants = []
        self.coverage_data = defaultdict(list)
        
    def load_variants(self):
        """Load variants from VCF file"""
        if not self.vcf_file or not os.path.exists(self.vcf_file):
            print("No VCF file provided or file doesn't exist")
            return
            
        vcf_reader = vcf.Reader(open(self.vcf_file, 'r'))
        
        for record in vcf_reader:
            variant = {
                'chrom': record.CHROM,
                'pos': record.POS,
                'ref': record.REF,
                'alt': str(record.ALT[0]) if record.ALT else '',
                'qual': record.QUAL,
                'filter': record.FILTER,
                'info': record.INFO,
                'type': self._classify_variant_type(record.REF, record.ALT[0] if record.ALT else '')
            }
            
            # Add genotype information if available
            if record.samples:
                sample = record.samples[0]
                variant['gt'] = sample['GT'] if 'GT' in sample.data else None
                variant['dp'] = sample['DP'] if 'DP' in sample.data else None
                variant['ad'] = sample['AD'] if 'AD' in sample.data else None
            
            self.variants.append(variant)
    
    def _classify_variant_type(self, ref, alt):
        """Classify variant type based on REF and ALT alleles"""
        if not alt:
            return 'unknown'
        
        alt = str(alt)
        if len(ref) == len(alt) == 1:
            return 'SNP'
        elif len(ref) > len(alt):
            return 'deletion'
        elif len(ref) < len(alt):
            return 'insertion'
        else:
            return 'complex'
    
    def analyze_coverage(self):
        """Analyze coverage from BAM file"""
        if not self.bam_file or not os.path.exists(self.bam_file):
            print("No BAM file provided or file doesn't exist")
            return
            
        bamfile = pysam.AlignmentFile(self.bam_file, "rb")
        
        for ref_name in self.references:
            ref_length = len(self.references[ref_name])
            
            # Get coverage for each position
            coverage = []
            for pos in range(ref_length):
                cov = bamfile.count_coverage(ref_name, pos, pos + 1)
                total_cov = sum(base_cov[0] for base_cov in cov)
                coverage.append(total_cov)
            
            self.coverage_data[ref_name] = coverage
        
        bamfile.close()
    
    def detect_structural_variants_from_alignment(self):
        """Detect structural variants from alignment data"""
        if not self.bam_file or not os.path.exists(self.bam_file):
            return
            
        bamfile = pysam.AlignmentFile(self.bam_file, "rb")
        
        for read in bamfile:
            if read.is_unmapped or read.is_secondary:
                continue
                
            # Analyze CIGAR string for structural variants
            self._analyze_cigar_for_sv(read)
        
        bamfile.close()
    
    def _analyze_cigar_for_sv(self, read):
        """Analyze CIGAR string to identify structural variants"""
        if not read.cigartuples:
            return
            
        ref_pos = read.reference_start
        query_pos = 0
        
        for op, length in read.cigartuples:
            if op == 1 and length >= 50:  # Large insertion
                sv = {
                    'type': 'insertion',
                    'chrom': read.reference_name,
                    'pos': ref_pos,
                    'length': length,
                    'sequence': read.query_sequence[query_pos:query_pos+length] if read.query_sequence else '',
                    'read_id': read.query_name,
                    'mapq': read.mapping_quality
                }
                self.structural_variants.append(sv)
                
            elif op == 2 and length >= 50:  # Large deletion
                sv = {
                    'type': 'deletion',
                    'chrom': read.reference_name,
                    'pos': ref_pos,
                    'length': length,
                    'sequence': '',
                    'read_id': read.query_name,
                    'mapq': read.mapping_quality
                }
                self.structural_variants.append(sv)
                
            # Update positions
            if op in [0, 2, 3]:  # M, D, N
                ref_pos += length
            if op in [0, 1, 4, 5]:  # M, I, S, H
                query_pos += length
    
    def calculate_mutation_density(self, window_size=100):
        """Calculate mutation density across reference sequences"""
        mutation_density = {}
        
        for ref_name in self.references:
            ref_length = len(self.references[ref_name])
            density = []
            
            # Count mutations in sliding windows
            for start in range(0, ref_length, window_size):
                end = min(start + window_size, ref_length)
                mutations_in_window = sum(1 for var in self.variants 
                                        if var['chrom'] == ref_name and start <= var['pos'] < end)
                density.append(mutations_in_window / window_size * 1000)  # Per kb
            
            mutation_density[ref_name] = {
                'positions': list(range(0, ref_length, window_size)),
                'density': density
            }
        
        return mutation_density
    
    def identify_hotspots(self, window_size=50, threshold=5):
        """Identify mutation hotspots"""
        hotspots = []
        
        for ref_name in self.references:
            ref_length = len(self.references[ref_name])
            
            for start in range(0, ref_length, window_size):
                end = min(start + window_size, ref_length)
                mutations_in_window = [var for var in self.variants 
                                     if var['chrom'] == ref_name and start <= var['pos'] < end]
                
                if len(mutations_in_window) >= threshold:
                    hotspots.append({
                        'chrom': ref_name,
                        'start': start,
                        'end': end,
                        'mutation_count': len(mutations_in_window),
                        'mutations': mutations_in_window
                    })
        
        return hotspots
    
    def annotate_variants(self):
        """Annotate variants with functional information"""
        annotated_variants = []
        
        for variant in self.variants:
            ref_seq = self.references[variant['chrom']].seq
            pos = variant['pos'] - 1  # Convert to 0-based
            
            # Basic annotation
            annotation = {
                'variant': variant,
                'context': str(ref_seq[max(0, pos-10):pos+10]),
                'gc_content': self._calculate_gc_content(ref_seq[max(0, pos-50):pos+50]),
            }
            
            # Predict effect (simplified)
            annotation['predicted_effect'] = self._predict_variant_effect(variant, ref_seq, pos)
            
            annotated_variants.append(annotation)
        
        return annotated_variants
    
    def _calculate_gc_content(self, sequence):
        """Calculate GC content of a sequence"""
        if not sequence:
            return 0
        return (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
    
    def _predict_variant_effect(self, variant, ref_seq, pos):
        """Predict the effect of a variant (simplified prediction)"""
        if variant['type'] == 'SNP':
            return 'substitution'
        elif variant['type'] == 'insertion':
            if len(variant['alt']) % 3 == 0:
                return 'in-frame_insertion'
            else:
                return 'frameshift_insertion'
        elif variant['type'] == 'deletion':
            if len(variant['ref']) % 3 == 0:
                return 'in-frame_deletion'
            else:
                return 'frameshift_deletion'
        else:
            return 'complex_variant'
    
    def generate_visualizations(self):
        """Generate comprehensive visualizations"""
        # Set up plotting style
        plt.style.use('seaborn-v0_8')
        fig_size = (15, 10)
        
        # 1. Variant type distribution
        if self.variants:
            self._plot_variant_distribution()
        
        # 2. Coverage plots
        if self.coverage_data:
            self._plot_coverage()
        
        # 3. Mutation density
        if self.variants:
            self._plot_mutation_density()
        
        # 4. Quality distribution
        if self.variants:
            self._plot_quality_distribution()
        
        # 5. Structural variant analysis
        if self.structural_variants:
            self._plot_structural_variants()
    
    def _plot_variant_distribution(self):
        """Plot distribution of variant types"""
        variant_types = [var['type'] for var in self.variants]
        type_counts = Counter(variant_types)
        
        plt.figure(figsize=(10, 6))
        plt.bar(type_counts.keys(), type_counts.values())
        plt.title('Distribution of Variant Types')
        plt.xlabel('Variant Type')
        plt.ylabel('Count')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'variant_type_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _plot_coverage(self):
        """Plot coverage across reference sequences"""
        n_refs = len(self.coverage_data)
        fig, axes = plt.subplots(n_refs, 1, figsize=(15, 4*n_refs))
        if n_refs == 1:
            axes = [axes]
        
        for i, (ref_name, coverage) in enumerate(self.coverage_data.items()):
            axes[i].plot(coverage)
            axes[i].set_title(f'Coverage: {ref_name}')
            axes[i].set_xlabel('Position (bp)')
            axes[i].set_ylabel('Coverage')
            axes[i].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'coverage_plots.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _plot_mutation_density(self):
        """Plot mutation density across sequences"""
        density_data = self.calculate_mutation_density()
        
        n_refs = len(density_data)
        fig, axes = plt.subplots(n_refs, 1, figsize=(15, 4*n_refs))
        if n_refs == 1:
            axes = [axes]
        
        for i, (ref_name, data) in enumerate(density_data.items()):
            axes[i].plot(data['positions'], data['density'])
            axes[i].set_title(f'Mutation Density: {ref_name}')
            axes[i].set_xlabel('Position (bp)')
            axes[i].set_ylabel('Mutations per kb')
            axes[i].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'mutation_density.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _plot_quality_distribution(self):
        """Plot quality score distribution of variants"""
        qualities = [var['qual'] for var in self.variants if var['qual'] is not None]
        
        if not qualities:
            return
        
        plt.figure(figsize=(10, 6))
        plt.hist(qualities, bins=50, alpha=0.7, edgecolor='black')
        plt.title('Variant Quality Score Distribution')
        plt.xlabel('Quality Score')
        plt.ylabel('Count')
        plt.axvline(np.mean(qualities), color='red', linestyle='--', label=f'Mean: {np.mean(qualities):.1f}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'quality_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _plot_structural_variants(self):
        """Plot structural variant analysis"""
        if not self.structural_variants:
            return
        
        sv_types = [sv['type'] for sv in self.structural_variants]
        sv_lengths = [sv['length'] for sv in self.structural_variants]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # SV type distribution
        type_counts = Counter(sv_types)
        ax1.bar(type_counts.keys(), type_counts.values())
        ax1.set_title('Structural Variant Types')
        ax1.set_xlabel('SV Type')
        ax1.set_ylabel('Count')
        
        # SV length distribution
        ax2.hist(sv_lengths, bins=30, alpha=0.7, edgecolor='black')
        ax2.set_title('Structural Variant Length Distribution')
        ax2.set_xlabel('Length (bp)')
        ax2.set_ylabel('Count')
        ax2.set_yscale('log')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'structural_variants.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def generate_summary_report(self):
        """Generate comprehensive summary report"""
        # Calculate statistics
        total_variants = len(self.variants)
        variant_types = Counter(var['type'] for var in self.variants)
        total_sv = len(self.structural_variants)
        sv_types = Counter(sv['type'] for sv in self.structural_variants)
        
        # Coverage statistics
        coverage_stats = {}
        for ref_name, coverage in self.coverage_data.items():
            if coverage:
                coverage_stats[ref_name] = {
                    'mean': np.mean(coverage),
                    'median': np.median(coverage),
                    'min': min(coverage),
                    'max': max(coverage),
                    'std': np.std(coverage)
                }
        
        # Generate report
        report = {
            'summary': {
                'total_variants': total_variants,
                'total_structural_variants': total_sv,
                'reference_sequences': len(self.references)
            },
            'variant_types': dict(variant_types),
            'structural_variant_types': dict(sv_types),
            'coverage_statistics': coverage_stats,
            'hotspots': self.identify_hotspots()
        }
        
        # Write JSON report
        with open(self.output_dir / 'analysis_summary.json', 'w') as f:
            json.dump(report, f, indent=2, default=str)
        
        # Write human-readable report
        with open(self.output_dir / 'analysis_report.txt', 'w') as f:
            f.write("PLASMID VARIANT ANALYSIS REPORT\n")
            f.write("=" * 40 + "\n\n")
            
            f.write("SUMMARY:\n")
            f.write(f"Total variants: {total_variants}\n")
            f.write(f"Total structural variants: {total_sv}\n")
            f.write(f"Reference sequences: {len(self.references)}\n\n")
            
            f.write("VARIANT TYPES:\n")
            for vtype, count in variant_types.items():
                f.write(f"  {vtype}: {count}\n")
            f.write("\n")
            
            if sv_types:
                f.write("STRUCTURAL VARIANT TYPES:\n")
                for svtype, count in sv_types.items():
                    f.write(f"  {svtype}: {count}\n")
                f.write("\n")
            
            if coverage_stats:
                f.write("COVERAGE STATISTICS:\n")
                for ref_name, stats in coverage_stats.items():
                    f.write(f"  {ref_name}:\n")
                    f.write(f"    Mean coverage: {stats['mean']:.1f}x\n")
                    f.write(f"    Median coverage: {stats['median']:.1f}x\n")
                    f.write(f"    Coverage range: {stats['min']}-{stats['max']}x\n")
        
        return report
    
    def run_complete_analysis(self):
        """Run complete variant analysis pipeline"""
        print("Starting variant analysis...")
        
        print("Loading variants...")
        self.load_variants()
        
        print("Analyzing coverage...")
        self.analyze_coverage()
        
        print("Detecting structural variants...")
        self.detect_structural_variants_from_alignment()
        
        print("Generating visualizations...")
        self.generate_visualizations()
        
        print("Generating summary report...")
        report = self.generate_summary_report()
        
        print(f"Analysis complete! Results saved to {self.output_dir}")
        return report


def main():
    parser = argparse.ArgumentParser(description="Advanced Plasmid Variant Analysis")
    parser.add_argument("--reference", "-r", required=True, help="Reference FASTA file")
    parser.add_argument("--vcf", "-v", help="VCF file with variants")
    parser.add_argument("--bam", "-b", help="BAM file with alignments")
    parser.add_argument("--output", "-o", default="variant_analysis", help="Output directory")
    parser.add_argument("--window-size", type=int, default=100, help="Window size for mutation density")
    parser.add_argument("--hotspot-threshold", type=int, default=5, help="Minimum mutations for hotspot")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.reference):
        print(f"Error: Reference file {args.reference} not found")
        sys.exit(1)
    
    if args.vcf and not os.path.exists(args.vcf):
        print(f"Warning: VCF file {args.vcf} not found")
    
    if args.bam and not os.path.exists(args.bam):
        print(f"Warning: BAM file {args.bam} not found")
    
    # Run analysis
    analyzer = PlasmidVariantAnalyzer(
        reference_file=args.reference,
        vcf_file=args.vcf,
        bam_file=args.bam,
        output_dir=args.output
    )
    
    analyzer.run_complete_analysis()
    print("Analysis completed successfully!")


if __name__ == "__main__":
    main()