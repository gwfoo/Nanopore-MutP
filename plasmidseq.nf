#!/usr/bin/env nextflow

/*
 * Plasmid Assembly and Variant Calling Pipeline
 * Oxford Nanopore MinION Data Processing
 */

nextflow.enable.dsl = 2

// Parameters
params.reads = "data/*.fastq.gz"
params.reference = "reference.fasta"
params.min_plasmid_size = 2000
params.max_plasmid_size = 20000
params.outdir = "results"
params.blast_db = "/path/to/blast/nt"
params.min_read_length = 1000
params.max_read_length = 50000
params.min_coverage = 10
params.overlap_threshold = 1000

// Tool versions and containers
params.flye_container = "staphb/flye:2.9.2"
params.minimap2_container = "staphb/minimap2:2.24"
params.samtools_container = "staphb/samtools:1.17"
params.medaka_container = "ontresearch/medaka:1.8.0"
params.blast_container = "ncbi/blast:2.14.0"
params.python_container = "python:3.9-slim"

workflow {
    // Input channels
    reads_ch = Channel.fromPath(params.reads, checkIfExists: true)
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true)
    
    // Main workflow
    QUALITY_FILTER(reads_ch)
    SIZE_FILTER(QUALITY_FILTER.out.filtered_reads, params.min_plasmid_size, params.max_plasmid_size)
    PLASMID_ASSEMBLY(SIZE_FILTER.out.plasmid_reads)
    POLISH_ASSEMBLY(PLASMID_ASSEMBLY.out.assembly, SIZE_FILTER.out.plasmid_reads)
    IDENTIFY_PLASMIDS(POLISH_ASSEMBLY.out.polished_assembly, params.min_plasmid_size, params.max_plasmid_size)
    ALIGN_TO_REFERENCE(IDENTIFY_PLASMIDS.out.plasmids, reference_ch)
    VARIANT_CALLING(ALIGN_TO_REFERENCE.out.alignments, reference_ch)
    STRUCTURAL_VARIANT_DETECTION(ALIGN_TO_REFERENCE.out.alignments, reference_ch)
    BLAST_UNKNOWN_REGIONS(STRUCTURAL_VARIANT_DETECTION.out.insertions)
    GENERATE_REPORT(
        VARIANT_CALLING.out.variants,
        STRUCTURAL_VARIANT_DETECTION.out.sv_results,
        BLAST_UNKNOWN_REGIONS.out.blast_results,
        IDENTIFY_PLASMIDS.out.plasmid_info
    )
}

process QUALITY_FILTER {
    container params.python_container
    publishDir "${params.outdir}/01_quality_filter", mode: 'copy'
    
    input:
    path reads
    
    output:
    path "filtered_reads.fastq.gz", emit: filtered_reads
    path "quality_stats.txt", emit: stats
    
    script:
    """
    #!/usr/bin/env python3
    import gzip
    from Bio import SeqIO
    import statistics
    
    def filter_reads(input_file, output_file, min_length=${params.min_read_length}, max_length=${params.max_read_length}, min_qual=7):
        filtered_count = 0
        total_count = 0
        lengths = []
        qualities = []
        
        with gzip.open(input_file, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
            for record in SeqIO.parse(infile, 'fastq'):
                total_count += 1
                read_length = len(record.seq)
                avg_qual = statistics.mean(record.letter_annotations['phred_quality'])
                
                if min_length <= read_length <= max_length and avg_qual >= min_qual:
                    SeqIO.write(record, outfile, 'fastq')
                    filtered_count += 1
                    lengths.append(read_length)
                    qualities.append(avg_qual)
        
        # Write statistics
        with open('quality_stats.txt', 'w') as stats:
            stats.write(f"Total reads: {total_count}\\n")
            stats.write(f"Filtered reads: {filtered_count}\\n")
            stats.write(f"Pass rate: {filtered_count/total_count*100:.2f}%\\n")
            if lengths:
                stats.write(f"Mean length: {statistics.mean(lengths):.0f}\\n")
                stats.write(f"Median length: {statistics.median(lengths):.0f}\\n")
                stats.write(f"Mean quality: {statistics.mean(qualities):.2f}\\n")
    
    filter_reads('${reads}', 'filtered_reads.fastq.gz')
    """
}

process SIZE_FILTER {
    container params.minimap2_container
    publishDir "${params.outdir}/02_size_filter", mode: 'copy'
    
    input:
    path reads
    val min_size
    val max_size
    
    output:
    path "plasmid_reads.fastq.gz", emit: plasmid_reads
    path "size_distribution.txt", emit: size_stats
    
    script:
    """
    # Use minimap2 to self-align reads and identify circular/plasmid-like sequences
    minimap2 -ax ava-ont ${reads} ${reads} > self_align.sam
    
    # Filter reads based on self-alignment patterns indicating circularity
    python3 << 'EOF'
import pysam
import gzip
from Bio import SeqIO
from collections import defaultdict

# Analyze self-alignments for circular patterns
circular_reads = set()
with pysam.AlignmentFile('self_align.sam', 'r') as samfile:
    for read in samfile:
        if not read.is_unmapped and read.reference_name == read.query_name:
            # Look for reads that align to themselves with wraparound
            if (read.query_alignment_start < 100 and 
                read.query_alignment_end > read.query_length - 100):
                circular_reads.add(read.query_name)

# Filter reads by size and circularity indicators
size_counts = defaultdict(int)
with gzip.open('${reads}', 'rt') as infile, gzip.open('plasmid_reads.fastq.gz', 'wt') as outfile:
    for record in SeqIO.parse(infile, 'fastq'):
        read_length = len(record.seq)
        size_counts[read_length//1000] += 1
        
        # Include reads in plasmid size range or showing circular patterns
        if (${min_size} <= read_length <= ${max_size} or 
            record.id in circular_reads):
            SeqIO.write(record, outfile, 'fastq')

# Write size distribution
with open('size_distribution.txt', 'w') as f:
    f.write("Size_kb\\tCount\\n")
    for size_kb, count in sorted(size_counts.items()):
        f.write(f"{size_kb}\\t{count}\\n")
EOF
    """
}

process PLASMID_ASSEMBLY {
    container params.flye_container
    publishDir "${params.outdir}/03_assembly", mode: 'copy'
    cpus 8
    memory '16 GB'
    
    input:
    path reads
    
    output:
    path "assembly/assembly.fasta", emit: assembly
    path "assembly/assembly_info.txt", emit: assembly_info
    path "flye_log.txt", emit: log
    
    script:
    """
    # Use Flye for de novo assembly optimized for plasmids
    flye --nano-raw ${reads} --out-dir assembly --threads ${task.cpus} \\
         --plasmids --meta --min-overlap ${params.overlap_threshold} > flye_log.txt 2>&1
    
    # Ensure output exists
    if [ ! -f assembly/assembly.fasta ]; then
        touch assembly/assembly.fasta
        echo "No assembly produced" > assembly/assembly_info.txt
    fi
    """
}

process POLISH_ASSEMBLY {
    container params.medaka_container
    publishDir "${params.outdir}/04_polishing", mode: 'copy'
    cpus 4
    memory '8 GB'
    
    input:
    path assembly
    path reads
    
    output:
    path "polished_assembly.fasta", emit: polished_assembly
    path "medaka_log.txt", emit: log
    
    script:
    """
    # Polish assembly with Medaka
    medaka_consensus -i ${reads} -d ${assembly} -o medaka_output \\
                     -t ${task.cpus} -m r941_min_high_g360 > medaka_log.txt 2>&1
    
    # Copy polished consensus
    if [ -f medaka_output/consensus.fasta ]; then
        cp medaka_output/consensus.fasta polished_assembly.fasta
    else
        cp ${assembly} polished_assembly.fasta
        echo "Polishing failed, using original assembly" >> medaka_log.txt
    fi
    """
}

process IDENTIFY_PLASMIDS {
    container params.python_container
    publishDir "${params.outdir}/05_plasmid_identification", mode: 'copy'
    
    input:
    path assembly
    val min_size
    val max_size
    
    output:
    path "plasmids.fasta", emit: plasmids
    path "plasmid_info.txt", emit: plasmid_info
    
    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import re
    
    def is_circular(sequence):
        \"\"\"Check for potential circularity by looking for overlap at ends\"\"\"
        seq_str = str(sequence).upper()
        end_len = min(100, len(seq_str) // 10)
        
        if end_len < 20:
            return False
            
        start_seq = seq_str[:end_len]
        end_seq = seq_str[-end_len:]
        
        # Check for overlap between start and end
        for i in range(20, end_len):
            if start_seq.endswith(end_seq[:i]) or end_seq.startswith(start_seq[-i:]):
                return True
        return False
    
    plasmids = []
    with open('${assembly}', 'r') as infile:
        for record in SeqIO.parse(infile, 'fasta'):
            length = len(record.seq)
            
            # Filter by size and check for plasmid characteristics
            if ${min_size} <= length <= ${max_size}:
                # Additional checks for plasmid-like sequences
                circular = is_circular(record.seq)
                gc_content = (record.seq.count('G') + record.seq.count('C')) / length * 100
                
                # Add plasmid if it meets criteria
                if circular or 'plasmid' in record.description.lower():
                    plasmids.append((record, length, circular, gc_content))
    
    # Write identified plasmids
    with open('plasmids.fasta', 'w') as outfile, open('plasmid_info.txt', 'w') as info:
        info.write("Plasmid_ID\\tLength\\tCircular\\tGC_Content\\n")
        
        for i, (record, length, circular, gc_content) in enumerate(plasmids):
            new_id = f"plasmid_{i+1}"
            record.id = new_id
            record.description = f"length={length} circular={circular}"
            SeqIO.write(record, outfile, 'fasta')
            info.write(f"{new_id}\\t{length}\\t{circular}\\t{gc_content:.2f}\\n")
    
    # Create empty files if no plasmids found
    if not plasmids:
        open('plasmids.fasta', 'w').close()
        with open('plasmid_info.txt', 'w') as info:
            info.write("Plasmid_ID\\tLength\\tCircular\\tGC_Content\\n")
            info.write("No plasmids identified\\n")
    """
}

process ALIGN_TO_REFERENCE {
    container params.minimap2_container
    publishDir "${params.outdir}/06_alignment", mode: 'copy'
    
    input:
    path plasmids
    path reference
    
    output:
    path "*.bam", emit: alignments
    path "*.bam.bai", emit: indices
    path "alignment_stats.txt", emit: stats
    
    script:
    """
    # Align each plasmid to reference
    minimap2 -ax asm5 ${reference} ${plasmids} | samtools sort -o plasmids_to_ref.bam
    samtools index plasmids_to_ref.bam
    
    # Generate alignment statistics
    samtools flagstat plasmids_to_ref.bam > alignment_stats.txt
    samtools coverage plasmids_to_ref.bam >> alignment_stats.txt
    """
}

process VARIANT_CALLING {
    container params.samtools_container
    publishDir "${params.outdir}/07_variants", mode: 'copy'
    
    input:
    path alignments
    path reference
    
    output:
    path "variants.vcf", emit: variants
    path "variant_summary.txt", emit: summary
    
    script:
    """
    # Call variants using bcftools
    samtools mpileup -uf ${reference} ${alignments} | \\
    bcftools call -mv -Ov -o raw_variants.vcf
    
    # Filter variants
    bcftools filter -s LOWQUAL -e 'QUAL<20 || DP<${params.min_coverage}' \\
                    raw_variants.vcf > variants.vcf
    
    # Generate summary
    echo "Variant Summary:" > variant_summary.txt
    echo "=================" >> variant_summary.txt
    bcftools stats variants.vcf >> variant_summary.txt
    """
}

process STRUCTURAL_VARIANT_DETECTION {
    container params.python_container
    publishDir "${params.outdir}/08_structural_variants", mode: 'copy'
    
    input:
    path alignments
    path reference
    
    output:
    path "structural_variants.txt", emit: sv_results
    path "insertions.fasta", emit: insertions
    path "sv_summary.txt", emit: summary
    
    script:
    """
    #!/usr/bin/env python3
    import pysam
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    import re
    
    def detect_structural_variants(bam_file, reference_file):
        # Load reference
        ref_dict = SeqIO.to_dict(SeqIO.parse(reference_file, 'fasta'))
        
        insertions = []
        deletions = []
        rearrangements = []
        
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            for read in bam:
                if read.is_unmapped or read.is_secondary:
                    continue
                
                # Parse CIGAR for structural variants
                cigar = read.cigartuples
                ref_pos = read.reference_start
                query_pos = 0
                
                for op, length in cigar:
                    if op == 1:  # Insertion
                        if length >= 50:  # Large insertion threshold
                            ins_seq = read.query_sequence[query_pos:query_pos+length]
                            insertions.append({
                                'contig': read.reference_name,
                                'position': ref_pos,
                                'length': length,
                                'sequence': ins_seq,
                                'read_id': read.query_name
                            })
                        query_pos += length
                    elif op == 2:  # Deletion
                        if length >= 50:  # Large deletion threshold
                            deletions.append({
                                'contig': read.reference_name,
                                'position': ref_pos,
                                'length': length,
                                'read_id': read.query_name
                            })
                        ref_pos += length
                    elif op == 0:  # Match
                        ref_pos += length
                        query_pos += length
                    elif op == 4:  # Soft clip
                        query_pos += length
        
        return insertions, deletions, rearrangements
    
    # Detect structural variants
    insertions, deletions, rearrangements = detect_structural_variants('${alignments}', '${reference}')
    
    # Write results
    with open('structural_variants.txt', 'w') as f:
        f.write("Type\\tContig\\tPosition\\tLength\\tRead_ID\\tSequence\\n")
        
        for ins in insertions:
            f.write(f"INS\\t{ins['contig']}\\t{ins['position']}\\t{ins['length']}\\t{ins['read_id']}\\t{ins['sequence'][:100]}\\n")
        
        for dele in deletions:
            f.write(f"DEL\\t{dele['contig']}\\t{dele['position']}\\t{dele['length']}\\t{dele['read_id']}\\tN/A\\n")
    
    # Write insertion sequences for BLAST
    with open('insertions.fasta', 'w') as f:
        for i, ins in enumerate(insertions):
            if len(ins['sequence']) >= 100:  # Only BLAST longer insertions
                record = SeqRecord(
                    Seq(ins['sequence']),
                    id=f"insertion_{i+1}",
                    description=f"pos:{ins['position']} len:{ins['length']}"
                )
                SeqIO.write(record, f, 'fasta')
    
    # Summary
    with open('sv_summary.txt', 'w') as f:
        f.write(f"Structural Variant Summary\\n")
        f.write(f"=========================\\n")
        f.write(f"Large insertions (≥50bp): {len(insertions)}\\n")
        f.write(f"Large deletions (≥50bp): {len(deletions)}\\n")
        f.write(f"Rearrangements: {len(rearrangements)}\\n")
    """
}

process BLAST_UNKNOWN_REGIONS {
    container params.blast_container
    publishDir "${params.outdir}/09_blast_analysis", mode: 'copy'
    
    input:
    path insertions
    
    output:
    path "blast_results.txt", emit: blast_results
    path "blast_summary.txt", emit: summary
    
    when:
    insertions.size() > 0
    
    script:
    """
    # BLAST insertions against nucleotide database
    if [ -s ${insertions} ]; then
        blastn -query ${insertions} -db ${params.blast_db} \\
               -out blast_results.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \\
               -max_target_seqs 5 -evalue 1e-5
        
        # Generate summary
        echo "BLAST Analysis Summary" > blast_summary.txt
        echo "=====================" >> blast_summary.txt
        echo "Query sequences: \$(grep -c '^>' ${insertions})" >> blast_summary.txt
        echo "Significant hits (E-value < 1e-5): \$(wc -l < blast_results.txt)" >> blast_summary.txt
        
        # Top hits summary
        echo "\\nTop hits:" >> blast_summary.txt
        head -10 blast_results.txt >> blast_summary.txt
    else
        echo "No insertions to BLAST" > blast_results.txt
        echo "No insertions found for BLAST analysis" > blast_summary.txt
    fi
    """
}

process GENERATE_REPORT {
    container params.python_container
    publishDir "${params.outdir}/10_final_report", mode: 'copy'
    
    input:
    path variants
    path sv_results
    path blast_results
    path plasmid_info
    
    output:
    path "final_report.html", emit: html_report
    path "summary_statistics.txt", emit: summary
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    from datetime import datetime
    import os
    
    def generate_html_report():
        html_content = f'''
    <!DOCTYPE html>
    <html>
    <head>
        <title>Plasmid Analysis Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            h1, h2 {{ color: #333; }}
            table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
            th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
            .summary {{ background-color: #f9f9f9; padding: 20px; border-radius: 5px; }}
        </style>
    </head>
    <body>
        <h1>Plasmid Assembly and Variant Analysis Report</h1>
        <p>Generated on: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
        
        <div class="summary">
            <h2>Summary</h2>
    '''
        
        # Add plasmid information
        try:
            with open('${plasmid_info}', 'r') as f:
                plasmid_lines = f.readlines()
                html_content += f"<p><strong>Identified Plasmids:</strong> {len(plasmid_lines)-1}</p>"
        except:
            html_content += "<p><strong>Identified Plasmids:</strong> 0</p>"
        
        # Add variant count
        try:
            with open('${variants}', 'r') as f:
                variant_count = sum(1 for line in f if not line.startswith('#'))
                html_content += f"<p><strong>Total Variants:</strong> {variant_count}</p>"
        except:
            html_content += "<p><strong>Total Variants:</strong> 0</p>"
        
        # Add structural variant count
        try:
            with open('${sv_results}', 'r') as f:
                sv_count = sum(1 for line in f) - 1  # Subtract header
                html_content += f"<p><strong>Structural Variants:</strong> {sv_count}</p>"
        except:
            html_content += "<p><strong>Structural Variants:</strong> 0</p>"
        
        html_content += '''
        </div>
        
        <h2>Detailed Results</h2>
        <p>Please refer to the individual output files for detailed analysis results:</p>
        <ul>
            <li>Variants: variants.vcf</li>
            <li>Structural Variants: structural_variants.txt</li>
            <li>BLAST Results: blast_results.txt</li>
            <li>Plasmid Information: plasmid_info.txt</li>
        </ul>
        
    </body>
    </html>
        '''
        
        with open('final_report.html', 'w') as f:
            f.write(html_content)
    
    def generate_summary_stats():
        with open('summary_statistics.txt', 'w') as f:
            f.write("Plasmid Analysis Summary Statistics\\n")
            f.write("===================================\\n\\n")
            
            # Plasmid stats
            try:
                with open('${plasmid_info}', 'r') as pf:
                    lines = pf.readlines()
                    f.write(f"Identified plasmids: {len(lines)-1}\\n")
            except:
                f.write("Identified plasmids: 0\\n")
            
            # Variant stats
            try:
                with open('${variants}', 'r') as vf:
                    variant_count = sum(1 for line in vf if not line.startswith('#'))
                    f.write(f"Total variants: {variant_count}\\n")
            except:
                f.write("Total variants: 0\\n")
            
            # SV stats
            try:
                with open('${sv_results}', 'r') as svf:
                    sv_count = sum(1 for line in svf) - 1
                    f.write(f"Structural variants: {sv_count}\\n")
            except:
                f.write("Structural variants: 0\\n")
    
    generate_html_report()
    generate_summary_stats()
    """
}

// Configuration
process {
    errorStrategy = 'retry'
    maxRetries = 2
    
    withLabel: 'high_memory' {
        memory = '32 GB'
        cpus = 8
    }
    
    withLabel: 'medium_memory' {
        memory = '16 GB'
        cpus = 4
    }
}

// Manifest
manifest {
    name = 'plasmid-assembly-pipeline'
    description = 'Oxford Nanopore plasmid assembly and variant calling pipeline'
    version = '1.0.0'
    mainScript = 'main.nf'
}