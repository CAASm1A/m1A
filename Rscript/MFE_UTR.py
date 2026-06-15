import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import os
import matplotlib.pyplot as plt
from scipy import stats
import argparse
import sys
from collections import defaultdict

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Analyze the impact of m1A modification on RNA secondary structure free energy',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze using full 5'UTR from GTF
  python m1a_mfe_analysis_full_utr.py --input m1a_sites.csv --fasta genome.fa --gtf annotation.gtf --output m1a_full_utr
  
  # Analyze using local window (original method)
  python m1a_mfe_analysis_full_utr.py --input m1a_sites.csv --fasta genome.fa --output m1a_window --window 100 --use-window
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input CSV file path (containing m1A site information)')
    
    parser.add_argument('-f', '--fasta', required=True,
                       help='Reference genome FASTA file path')
    
    parser.add_argument('-g', '--gtf', 
                       help='GTF annotation file path (required for full 5\'UTR analysis)')
    
    parser.add_argument('-o', '--output', default='m1a_mfe',
                       help='Output file prefix (default: m1a_mfe)')
    
    parser.add_argument('-w', '--window', type=int, default=100,
                       help='Window size for local analysis (default: 100)')
    
    parser.add_argument('--use-window', action='store_true',
                       help='Use local window analysis instead of full 5\'UTR')
    
    parser.add_argument('--max-utr-length', type=int, default=1000,
                       help='Maximum 5\'UTR length to analyze (default: 1000 bp)')
    
    parser.add_argument('--filter-utr', action='store_true', default=True,
                       help='Analyze only 5\'UTR regions (default: True)')
    
    parser.add_argument('--utr-column', type=int, default=5,
                       help='UTR information column index (0-based, default: 5)')
    
    parser.add_argument('--chrom-column', type=int, default=0,
                       help='Chromosome information column index (0-based, default: 0)')
    
    parser.add_argument('--pos-column', type=int, default=1,
                       help='Position information column index (1-based, default: 1)')
    
    parser.add_argument('--strand-column', type=int, default=10,
                       help='Strand information column index (default: 10)')
    
    parser.add_argument('--gene-column', type=int, default=11,
                       help='Gene name column index (default: 11)')
    
    return parser.parse_args()


def parse_gtf_attribute(attr_string, key):
    """Extract attribute value from GTF attribute string"""
    for item in attr_string.split(';'):
        item = item.strip()
        if item.startswith(key):
            return item.split('"')[1]
    return None


def load_5utr_from_gtf(gtf_file):
    """
    Load 5'UTR regions from GTF file
    Returns: dict of {gene_id: [(chrom, start, end, strand, transcript_id), ...]}
    """
    print(f"Loading 5'UTR regions from {gtf_file}...")
    
    utr_dict = defaultdict(list)
    transcript_utrs = defaultdict(list)  # {transcript_id: [(start, end), ...]}
    
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                chrom = fields[0]
                feature = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                attributes = fields[8]
                
                # Only process five_prime_utr features
                if feature != 'five_prime_utr':
                    continue
                
                gene_id = parse_gtf_attribute(attributes, 'gene_id')
                transcript_id = parse_gtf_attribute(attributes, 'transcript_id')
                
                if gene_id and transcript_id:
                    # Store by transcript first to merge exons
                    key = (chrom, strand, transcript_id, gene_id)
                    transcript_utrs[key].append((start, end))
        
        # Merge overlapping/adjacent UTR regions for each transcript
        for (chrom, strand, transcript_id, gene_id), regions in transcript_utrs.items():
            # Sort regions by start position
            regions.sort()
            
            # Merge overlapping regions
            merged = []
            current_start, current_end = regions[0]
            
            for start, end in regions[1:]:
                if start <= current_end + 1:  # Overlapping or adjacent
                    current_end = max(current_end, end)
                else:
                    merged.append((current_start, current_end))
                    current_start, current_end = start, end
            merged.append((current_start, current_end))
            
            # Store merged regions
            for start, end in merged:
                utr_dict[gene_id].append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'transcript_id': transcript_id,
                    'length': end - start + 1
                })
        
        print(f"  Loaded 5'UTR information for {len(utr_dict)} genes")
        total_utrs = sum(len(v) for v in utr_dict.values())
        print(f"  Total 5'UTR regions: {total_utrs}")
        
        return dict(utr_dict)
        
    except Exception as e:
        print(f"Error reading GTF file: {e}")
        sys.exit(1)


def extract_full_5utr_sequence(genome_dict, utr_info):
    """
    Extract full 5'UTR sequence from genome
    """
    chrom = utr_info['chrom']
    start = utr_info['start']
    end = utr_info['end']
    strand = utr_info['strand']
    
    if chrom not in genome_dict:
        return None
    
    # Extract sequence (GTF is 1-based, Python is 0-based)
    seq_obj = genome_dict[chrom].seq[start-1:end]
    
    # Handle negative strand
    if strand == '-':
        seq_str = str(seq_obj.reverse_complement())
    else:
        seq_str = str(seq_obj)
    
    # Convert to RNA (T -> U)
    seq_str = seq_str.replace('T', 'U')
    
    return seq_str


def find_m1a_position_in_utr(m1a_pos, utr_info):
    """
    Find the position of m1A within the 5'UTR sequence
    Returns: position in UTR (1-based), or None if not in this UTR
    """
    start = utr_info['start']
    end = utr_info['end']
    strand = utr_info['strand']
    
    # Check if m1A is within this UTR region
    if not (start <= m1a_pos <= end):
        return None
    
    if strand == '+':
        # For + strand, position from start
        pos_in_utr = m1a_pos - start + 1
    else:
        # For - strand, position from end (reverse)
        pos_in_utr = end - m1a_pos + 1
    
    return pos_in_utr


def calculate_mfe(sequence, constraint_pos=None):
    """
    Calculate MFE using RNAfold
    
    Parameters:
    sequence: RNA sequence
    constraint_pos: position to constrain (1-based), None for no constraint
    
    Returns: MFE value
    """
    try:
        if constraint_pos is None:
            # Native folding
            proc = subprocess.Popen(['RNAfold', '--noPS'], 
                                   stdin=subprocess.PIPE, 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE,
                                   text=True)
            out, _ = proc.communicate(input=sequence)
        else:
            # Constrained folding
            constraint = '.' * (constraint_pos - 1) + 'x' + '.' * (len(sequence) - constraint_pos)
            input_str = f"{sequence}\n{constraint}"
            
            proc = subprocess.Popen(['RNAfold', '-C', '--noPS'], 
                                   stdin=subprocess.PIPE, 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE,
                                   text=True)
            out, _ = proc.communicate(input=input_str)
        
        # Parse MFE from output
        mfe = float(out.split('(')[-1].split(')')[0])
        return mfe
        
    except Exception as e:
        print(f"Error in RNAfold: {e}")
        return None


def main():
    # --- Parse command line arguments ---
    args = parse_arguments()
    
    # Check required arguments
    if not args.use_window and not args.gtf:
        print("Error: --gtf is required for full 5'UTR analysis")
        print("Use --use-window flag for local window analysis without GTF")
        sys.exit(1)
    
    # Create output directory
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # --- Output parameters ---
    print("=" * 60)
    print("m1A MFE Analysis Tool - Full 5'UTR Version")
    print("=" * 60)
    print(f"Input file: {args.input}")
    print(f"Reference genome: {args.fasta}")
    if not args.use_window:
        print(f"GTF file: {args.gtf}")
        print(f"Analysis mode: Full 5'UTR")
        print(f"Max 5'UTR length: {args.max_utr_length} bp")
    else:
        print(f"Analysis mode: Local window (±{args.window} bp)")
    print(f"Output prefix: {args.output}")
    print("-" * 60)
    
    # --- Load genome ---
    print(f"Loading genome...")
    if not os.path.exists(args.fasta):
        print(f"Error: Genome file not found: {args.fasta}")
        sys.exit(1)
    
    try:
        genome_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))
        print(f"  Loaded {len(genome_dict)} chromosome(s)")
    except Exception as e:
        print(f"Error loading genome: {e}")
        sys.exit(1)
    
    # --- Load GTF if needed ---
    utr_dict = None
    if not args.use_window:
        utr_dict = load_5utr_from_gtf(args.gtf)
    
    # --- Read m1A sites ---
    print(f"Reading m1A sites...")
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)
    
    try:
        df = pd.read_csv(args.input, header=None)
        print(f"  Read {len(df)} rows")
    except Exception as e:
        print(f"Error reading CSV: {e}")
        sys.exit(1)
    
    # Filter 5'UTR
    if args.filter_utr:
        original_count = len(df)
        df = df[df[args.utr_column] == "5' UTR"]
        print(f"  Filtered to {len(df)} 5'UTR sites")
    
    # --- Process sites ---
    results = []
    processed = 0
    skipped = 0
    no_utr_found = 0
    utr_too_long = 0
    m1a_not_in_utr = 0
    
    print(f"Processing {len(df)} sites...")
    
    for index, row in df.iterrows():
        try:
            chrom = str(row[args.chrom_column])
            pos = int(row[args.pos_column])  # 1-based
            strand_code = row[args.strand_column]  # 1 for +, 2 for -
            gene_name = row[args.gene_column]
            
            strand = '+' if strand_code == 1 else '-'
            
            if args.use_window:
                # --- Local window analysis (original method) ---
                window_size = args.window
                start = pos - window_size - 1  # 0-based
                end = pos + window_size
                
                if chrom not in genome_dict:
                    skipped += 1
                    continue
                
                if start < 0 or end > len(genome_dict[chrom]):
                    skipped += 1
                    continue
                
                seq_obj = genome_dict[chrom].seq[start:end]
                
                if strand == '-':
                    seq_str = str(seq_obj.reverse_complement())
                else:
                    seq_str = str(seq_obj)
                
                seq_str = seq_str.replace('T', 'U')
                
                # m1A is at the center
                m1a_pos_in_seq = window_size + 1
                
            else:
                # --- Full 5'UTR analysis ---
                if gene_name not in utr_dict:
                    no_utr_found += 1
                    skipped += 1
                    continue
                
                # Find which UTR region contains this m1A
                found_utr = None
                m1a_pos_in_seq = None
                
                for utr_info in utr_dict[gene_name]:
                    if utr_info['chrom'] != chrom or utr_info['strand'] != strand:
                        continue
                    
                    pos_in_utr = find_m1a_position_in_utr(pos, utr_info)
                    if pos_in_utr is not None:
                        found_utr = utr_info
                        m1a_pos_in_seq = pos_in_utr
                        break
                
                if found_utr is None:
                    m1a_not_in_utr += 1
                    skipped += 1
                    continue
                
                # Check UTR length
                if found_utr['length'] > args.max_utr_length:
                    utr_too_long += 1
                    skipped += 1
                    continue
                
                # Extract full 5'UTR sequence
                seq_str = extract_full_5utr_sequence(genome_dict, found_utr)
                if seq_str is None:
                    skipped += 1
                    continue
            
            # --- Calculate MFE ---
            mfe_native = calculate_mfe(seq_str, constraint_pos=None)
            mfe_modified = calculate_mfe(seq_str, constraint_pos=m1a_pos_in_seq)
            
            if mfe_native is None or mfe_modified is None:
                skipped += 1
                continue
            
            result = {
                'Chr': chrom,
                'Position': pos,
                'Strand': strand,
                'Gene': gene_name,
                'Sequence_Length': len(seq_str),
                'm1A_Position_in_Seq': m1a_pos_in_seq,
                'Sequence': seq_str if len(seq_str) <= 200 else seq_str[:100] + '...' + seq_str[-100:],
                'MFE_Native': mfe_native,
                'MFE_Modified': mfe_modified,
                'Delta_MFE': mfe_modified - mfe_native
            }
            
            if not args.use_window:
                result['UTR_Start'] = found_utr['start']
                result['UTR_End'] = found_utr['end']
                result['Transcript_ID'] = found_utr['transcript_id']
            
            results.append(result)
            processed += 1
            
            if processed % 100 == 0:
                print(f"  Processed {processed} sites...")
                
        except Exception as e:
            print(f"  Error processing row {index}: {e}")
            skipped += 1
            continue
    
    print(f"\nProcessing complete:")
    print(f"  Successfully processed: {processed}")
    print(f"  Skipped: {skipped}")
    if not args.use_window:
        print(f"    - No UTR found: {no_utr_found}")
        print(f"    - UTR too long: {utr_too_long}")
        print(f"    - m1A not in UTR: {m1a_not_in_utr}")
    
    if len(results) == 0:
        print("\nError: No valid results generated")
        sys.exit(1)
    
    # --- Save results ---
    res_df = pd.DataFrame(results)
    output_csv = f"{args.output}_results.csv"
    res_df.to_csv(output_csv, index=False)
    print(f"\nResults saved: {output_csv}")
    
    # --- Statistics ---
    summary_file = f"{args.output}_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("m1A MFE Analysis Results Summary\n")
        f.write("=" * 60 + "\n")
        f.write(f"Input file: {args.input}\n")
        f.write(f"Reference genome: {args.fasta}\n")
        if not args.use_window:
            f.write(f"GTF file: {args.gtf}\n")
            f.write(f"Analysis mode: Full 5'UTR\n")
        else:
            f.write(f"Analysis mode: Local window (±{args.window} bp)\n")
        f.write(f"Total input sites: {len(df)}\n")
        f.write(f"Successfully processed: {processed}\n")
        f.write(f"Skipped: {skipped}\n")
        f.write("-" * 60 + "\n")
        f.write(f"Sequence length: {res_df['Sequence_Length'].mean():.1f} ± {res_df['Sequence_Length'].std():.1f} bp\n")
        f.write(f"Mean MFE (Native): {res_df['MFE_Native'].mean():.3f} ± {res_df['MFE_Native'].std():.3f} kcal/mol\n")
        f.write(f"Mean MFE (Modified): {res_df['MFE_Modified'].mean():.3f} ± {res_df['MFE_Modified'].std():.3f} kcal/mol\n")
        f.write(f"Mean ΔMFE: {res_df['Delta_MFE'].mean():.3f} ± {res_df['Delta_MFE'].std():.3f} kcal/mol\n")
        f.write(f"Median ΔMFE: {res_df['Delta_MFE'].median():.3f} kcal/mol\n")
        f.write("-" * 60 + "\n")
    
    # Statistical test
    if len(res_df) > 1:
        stat, pval = stats.ttest_rel(res_df['MFE_Native'], res_df['MFE_Modified'])
        
        with open(summary_file, 'a') as f:
            f.write(f"Paired t-test:\n")
            f.write(f"  t-statistic: {stat:.3f}\n")
            f.write(f"  p-value: {pval:.2e}\n")
            f.write("=" * 60 + "\n")
        
        # --- Plotting ---
        fig, axes = plt.subplots(2, 2, figsize=(14, 11))
        
        mean_native = res_df['MFE_Native'].mean()
        mean_modified = res_df['MFE_Modified'].mean()
        mean_delta = res_df['Delta_MFE'].mean()
        
        # Subplot 1: Box plot
        positions = [1, 2]
        boxes = axes[0, 0].boxplot([res_df['MFE_Native'], res_df['MFE_Modified']], 
                                   labels=['Native\n(No constraint)', 'With m1A\n(Constraint)'],
                                   positions=positions,
                                   patch_artist=True)
        
        colors = ['lightblue', 'lightcoral']
        for patch, color in zip(boxes['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        axes[0, 0].scatter([positions[0]], [mean_native], color='darkblue', s=100, 
                          marker='D', zorder=3, label=f'Mean: {mean_native:.1f}')
        axes[0, 0].scatter([positions[1]], [mean_modified], color='darkred', s=100, 
                          marker='D', zorder=3, label=f'Mean: {mean_modified:.1f}')
        axes[0, 0].plot(positions, [mean_native, mean_modified], 'k--', alpha=0.5, linewidth=1)
        axes[0, 0].set_ylabel('Minimum Free Energy (kcal/mol)')
        axes[0, 0].set_title('RNAfold MFE Comparison')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Subplot 2: ΔMFE distribution
        axes[0, 1].hist(res_df['Delta_MFE'], bins=30, edgecolor='black', alpha=0.7, color='lightgreen')
        axes[0, 1].axvline(mean_delta, color='red', linestyle='--', linewidth=2, 
                          label=f'Mean ΔMFE: {mean_delta:.2f}')
        axes[0, 1].axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
        axes[0, 1].set_xlabel('ΔMFE (Modified - Native, kcal/mol)')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('Distribution of ΔMFE Values')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # Subplot 3: Scatter plot
        axes[1, 0].scatter(res_df['MFE_Native'], res_df['MFE_Modified'], 
                          alpha=0.5, s=20, color='purple')
        
        min_val = min(res_df['MFE_Native'].min(), res_df['MFE_Modified'].min())
        max_val = max(res_df['MFE_Native'].max(), res_df['MFE_Modified'].max())
        axes[1, 0].plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, 
                       label='y = x')
        axes[1, 0].set_xlabel('MFE Native (kcal/mol)')
        axes[1, 0].set_ylabel('MFE Modified (kcal/mol)')
        axes[1, 0].set_title('MFE Native vs Modified')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # Subplot 4: Text summary
        axes[1, 1].axis('off')
        analysis_mode = "Full 5'UTR" if not args.use_window else f"Window ±{args.window}bp"
        summary_text = (
            f"Analysis Summary:\n\n"
            f"Mode: {analysis_mode}\n"
            f"Total sites: {len(res_df)}\n"
            f"Avg sequence length: {res_df['Sequence_Length'].mean():.0f} bp\n\n"
            f"Mean MFE (Native):\n  {mean_native:.2f} kcal/mol\n"
            f"Mean MFE (Modified):\n  {mean_modified:.2f} kcal/mol\n"
            f"Mean ΔMFE:\n  {mean_delta:.2f} ± {res_df['Delta_MFE'].std():.2f}\n"
            f"Median ΔMFE: {res_df['Delta_MFE'].median():.2f}\n\n"
            f"p-value: {pval:.2e}\n\n"
            f"ΔMFE > 0: Less stable\n"
            f"ΔMFE < 0: More stable"
        )
        axes[1, 1].text(0.1, 0.5, summary_text, fontsize=11, 
                       verticalalignment='center', fontfamily='monospace',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
        
        mode_title = "Full 5'UTR Analysis" if not args.use_window else f"Local Window Analysis (±{args.window}bp)"
        plt.suptitle(f'Impact of m1A on RNA Structure - {mode_title} (n={len(res_df)})', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        plot_file = f"{args.output}_plots.pdf"
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Summary saved: {summary_file}")
        print(f"Plot saved: {plot_file}")
        print(f"\nKey Results:")
        print(f"  Mean ΔMFE: {mean_delta:.3f} kcal/mol")
        print(f"  P-value: {pval:.2e}")
        print(f"  Average sequence length: {res_df['Sequence_Length'].mean():.1f} bp")
    
    print("=" * 60)
    print("Analysis completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()