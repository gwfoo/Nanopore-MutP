#!/bin/bash

# Plasmid Assembly Pipeline Setup and Execution Script
# This script sets up the environment and runs the plasmid analysis pipeline

set -e  # Exit on any error

# Configuration
PIPELINE_NAME="plasmid-assembly-pipeline"
CONDA_ENV_NAME="plasmid_pipeline"
PIPELINE_VERSION="1.0.0"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_header() {
    echo -e "${BLUE}$1${NC}"
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to install conda/mamba if not present
install_conda() {
    if ! command_exists conda; then
        print_status "Installing Miniconda..."
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
        bash miniconda.sh -b -p $HOME/miniconda
        export PATH="$HOME/miniconda/bin:$PATH"
        source $HOME/miniconda/etc/profile.d/conda.sh
        rm miniconda.sh
    fi
    
    # Install mamba for faster package management
    if ! command_exists mamba; then
        print_status "Installing Mamba..."
        conda install -y -c conda-forge mamba
    fi
}

# Function to create conda environment
create_environment() {
    print_status "Creating conda environment: $CONDA_ENV_NAME"
    
    cat > environment.yml << EOF
name: $CONDA_ENV_NAME
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.9
  - nextflow>=21.04
  - biopython
  - pysam
  - pandas
  - numpy
  - matplotlib
  - seaborn
  - pyvcf
  - samtools
  - bcftools
  - minimap2
  - blast
  - flye
  - medaka
  - pip
  - pip:
    - nextflow
EOF

    if conda env list | grep -q "$CONDA_ENV_NAME"; then
        print_warning "Environment $CONDA_ENV_NAME already exists. Updating..."
        mamba env update -n $CONDA_ENV_NAME -f environment.yml
    else
        mamba env create -f environment.yml
    fi
    
    rm environment.yml
}

# Function to setup directory structure
setup_directories() {
    print_status "Setting up directory structure..."
    
    mkdir -p {data,results,test_data,scripts,config}
    mkdir -p results/{01_quality_filter,02_size_filter,03_assembly,04_polishing,05_plasmid_identification,06_alignment,07_variants,08_structural_variants,09_blast_analysis,10_final_report}
    
    print_status "Directory structure created:"
    tree -d . 2>/dev/null || ls -la
}

# Function to download test data
download_test_data() {
    print_status "Setting up test data..."
    
    # Create synthetic test data if real data not available
    mkdir -p test_data
    
    # Generate a synthetic reference plasmid
    python3 << 'EOF'
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random

# Generate synthetic plasmid reference
def generate_plasmid_sequence(length=5000):
    bases = ['A', 'T', 'G', 'C']
    sequence = ''.join(random.choices(bases, k=length))
    return sequence

# Create reference plasmid
ref_seq = generate_plasmid_sequence(5000)
ref_record = SeqRecord(Seq(ref_seq), id="test_plasmid", description="Synthetic test plasmid")

with open("test_data/reference.fasta", "w") as f:
    SeqIO.write(ref_record, f, "fasta")

print("Synthetic reference plasmid created: test_data/reference.fasta")
EOF

    # Generate synthetic Oxford Nanopore reads
    python3 << 'EOF'
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random
import gzip

def generate_ont_reads(reference_seq, num_reads=1000, error_rate=0.10):
    reads = []
    ref_len = len(reference_seq)
    
    for i in range(num_reads):
        # Random read length (typical for ONT)
        read_length = random.randint(500, 8000)
        
        # Random start position (allow wraparound for circular plasmid)
        start_pos = random.randint(0, ref_len - 1)
        
        # Extract sequence (with potential wraparound)
        if start_pos + read_length <= ref_len:
            read_seq = reference_seq[start_pos:start_pos + read_length]
        else:
            # Wraparound for circular plasmid
            part1 = reference_seq[start_pos:]
            part2 = reference_seq[:read_length - len(part1)]
            read_seq = part1 + part2
        
        # Add sequencing errors
        read_seq = introduce_errors(read_seq, error_rate)
        
        # Generate quality scores (typical for ONT)
        qualities = [random.randint(5, 25) for _ in range(len(read_seq))]
        
        read_record = SeqRecord(
            Seq(read_seq),
            id=f"read_{i+1}",
            description=f"synthetic ONT read {i+1}"
        )
        read_record.letter_annotations["phred_quality"] = qualities
        reads.append(read_record)
    
    return reads

def introduce_errors(sequence, error_rate):
    seq_list = list(sequence)
    num_errors = int(len(sequence) * error_rate)
    
    for _ in range(num_errors):
        pos = random.randint(0, len(seq_list) - 1)
        error_type = random.choice(['substitution', 'insertion', 'deletion'])
        
        if error_type == 'substitution':
            seq_list[pos] = random.choice(['A', 'T', 'G', 'C'])
        elif error_type == 'insertion':
            seq_list.insert(pos, random.choice(['A', 'T', 'G', 'C']))
        elif error_type == 'deletion' and len(seq_list) > 1:
            seq_list.pop(pos)
    
    return ''.join(seq_list)

# Load reference and generate reads
with open("test_data/reference.fasta", "r") as f:
    ref_record = SeqIO.read(f, "fasta")

reads = generate_ont_reads(str(ref_record.seq), num_reads=500)

# Write reads to compressed FASTQ
with gzip.open("test_data/ont_reads.fastq.gz", "wt") as f:
    SeqIO.write(reads, f, "fastq")

print("Synthetic ONT reads created: test_data/ont_reads.fastq.gz")
EOF

    print_status "Test data created in test_data/ directory"
}

# Function to check system requirements
check_requirements() {
    print_header "Checking system requirements..."
    
    # Check available memory
    if command_exists free; then
        total_mem=$(free -g | awk '/^Mem:/{print $2}')
        if [ "$total_mem" -lt 8 ]; then
            print_warning "System has less than 8GB RAM. Pipeline may run slowly."
        else
            print_status "Memory: ${total_mem}GB available"
        fi
    fi
    
    # Check available disk space
    if command_exists df; then
        available_space=$(df -h . | awk 'NR==2{print $4}')
        print_status "Available disk space: $available_space"
    fi
    
    # Check for Docker/Singularity
    if command_exists docker; then
        print_status "Docker found: $(docker --version)"
    elif command_exists singularity; then
        print_status "Singularity found: $(singularity --version)"
    else
        print_warning "Neither Docker nor Singularity found. Installing containers may be slower."
    fi
}

# Function to validate input files
validate_inputs() {
    local reads_file=$1
    local reference_file=$2
    
    print_header "Validating input files..."
    
    if [ ! -z "$reads_file" ]; then
        if [ -f "$reads_file" ]; then
            print_status "Reads file found: $reads_file"
            # Check if file is compressed
            if [[ $reads_file == *.gz ]]; then
                read_count=$(zcat "$reads_file" | grep -c "^@" || echo "0")
            else
                read_count=$(grep -c "^@" "$reads_file" || echo "0")
            fi
            print_status "Estimated read count: $read_count"
        else
            print_error "Reads file not found: $reads_file"
            return 1
        fi
    fi
    
    if [ ! -z "$reference_file" ]; then
        if [ -f "$reference_file" ]; then
            print_status "Reference file found: $reference_file"
            seq_count=$(grep -c "^>" "$reference_file" || echo "0")
            print_status "Reference sequences: $seq_count"
        else
            print_error "Reference file not found: $reference_file"
            return 1
        fi
    fi
    
    return 0
}

# Function to run the pipeline
run_pipeline() {
    local reads_pattern=$1
    local reference_file=$2
    local output_dir=$3
    local profile=${4:-"docker"}
    
    print_header "Running Plasmid Assembly Pipeline..."
    
    # Activate conda environment
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate $CONDA_ENV_NAME
    
    # Check if Nextflow is available
    if ! command_exists nextflow; then
        print_error "Nextflow not found in environment"
        return 1
    fi
    
    # Create nextflow command
    local cmd="nextflow run main.nf"
    cmd="$cmd --reads '$reads_pattern'"
    cmd="$cmd --reference '$reference_file'"
    cmd="$cmd --outdir '$output_dir'"
    cmd="$cmd -profile $profile"
    cmd="$cmd -resume"
    
    print_status "Executing: $cmd"
    
    # Run pipeline
    eval $cmd
    
    if [ $? -eq 0 ]; then
        print_status "Pipeline completed successfully!"
        print_status "Results available in: $output_dir"
    else
        print_error "Pipeline failed. Check logs for details."
        return 1
    fi
}

# Function to run post-pipeline analysis
run_post_analysis() {
    local results_dir=$1
    local reference_file=$2
    
    print_header "Running post-pipeline variant analysis..."
    
    # Find generated files
    local vcf_file=$(find "$results_dir" -name "*.vcf" | head -1)
    local bam_file=$(find "$results_dir" -name "*.bam" | head -1)
    
    if [ -f "scripts/variant_analyzer.py" ]; then
        python scripts/variant_analyzer.py \
            --reference "$reference_file" \
            --vcf "$vcf_file" \
            --bam "$bam_file" \
            --output "$results_dir/variant_analysis"
        
        if [ $? -eq 0 ]; then
            print_status "Post-analysis completed successfully!"
            print_status "Detailed analysis available in: $results_dir/variant_analysis"
        else
            print_warning "Post-analysis failed, but main pipeline results are still available"
        fi
    else
        print_warning "Variant analyzer script not found, skipping detailed analysis"
    fi
}

# Function to generate final report
generate_final_report() {
    local results_dir=$1
    
    print_header "Generating final summary report..."
    
    cat > "$results_dir/PIPELINE_SUMMARY.md" << EOF
# Plasmid Assembly Pipeline Results

## Pipeline Information
- Pipeline: $PIPELINE_NAME
- Version: $PIPELINE_VERSION
- Run Date: $(date)
- Results Directory: $results_dir

## Output Structure
\`\`\`
$results_dir/
├── 01_quality_filter/     # Quality filtered reads
├── 02_size_filter/        # Size-filtered reads
├── 03_assembly/           # De novo assembly results
├── 04_polishing/          # Polished assembly
├── 05_plasmid_identification/  # Identified plasmids
├── 06_alignment/          # Reference alignment
├── 07_variants/           # Variant calling results
├── 08_structural_variants/     # Structural variants
├── 09_blast_analysis/     # BLAST results for unknown regions
├── 10_final_report/       # Final HTML report
└── variant_analysis/      # Detailed variant analysis
\`\`\`

## Key Results Files
- **Assembled Plasmids**: \`05_plasmid_identification/plasmids.fasta\`
- **Variants**: \`07_variants/variants.vcf\`
- **Structural Variants**: \`08_structural_variants/structural_variants.txt\`
- **BLAST Results**: \`09_blast_analysis/blast_results.txt\`
- **Final Report**: \`10_final_report/final_report.html\`

## Next Steps
1. Review the HTML report: \`10_final_report/final_report.html\`
2. Examine variant calls: \`07_variants/variants.vcf\`
3. Check structural variants: \`08_structural_variants/structural_variants.txt\`
4. Review detailed analysis: \`variant_analysis/analysis_report.txt\`

## Quality Metrics
- Check assembly statistics in \`03_assembly/assembly_info.txt\`
- Review alignment statistics in \`06_alignment/alignment_stats.txt\`
- Examine coverage plots in \`variant_analysis/coverage_plots.png\`

## Troubleshooting
If you encounter issues:
1. Check the Nextflow execution report: \`pipeline_report.html\`
2. Review process logs in the \`work/\` directory
3. Validate input file formats and paths
4. Ensure sufficient system resources (RAM, disk space)

EOF

    print_status "Pipeline summary created: $results_dir/PIPELINE_SUMMARY.md"
}

# Function to clean up temporary files
cleanup() {
    print_status "Cleaning up temporary files..."
    
    # Remove work directory if requested
    if [ "$CLEANUP_WORK" = "true" ]; then
        rm -rf work/
        print_status "Removed work directory"
    fi
    
    # Remove .nextflow directory if requested
    if [ "$CLEANUP_NEXTFLOW" = "true" ]; then
        rm -rf .nextflow*
        print_status "Removed .nextflow directories"
    fi
}

# Main execution function
main() {
    print_header "=== Plasmid Assembly Pipeline Setup ==="
    
    # Parse command line arguments
    SETUP_ONLY=false
    TEST_RUN=false
    READS_PATTERN=""
    REFERENCE_FILE=""
    OUTPUT_DIR="results"
    PROFILE="docker"
    CLEANUP_WORK=false
    CLEANUP_NEXTFLOW=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            --setup-only)
                SETUP_ONLY=true
                shift
                ;;
            --test)
                TEST_RUN=true
                shift
                ;;
            --reads)
                READS_PATTERN="$2"
                shift 2
                ;;
            --reference)
                REFERENCE_FILE="$2"
                shift 2
                ;;
            --output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --profile)
                PROFILE="$2"
                shift 2
                ;;
            --cleanup-work)
                CLEANUP_WORK=true
                shift
                ;;
            --cleanup-all)
                CLEANUP_WORK=true
                CLEANUP_NEXTFLOW=true
                shift
                ;;
            --help)
                cat << EOF
Usage: $0 [OPTIONS]

OPTIONS:
    --setup-only        Only setup environment, don't run pipeline
    --test             Run with test data
    --reads PATTERN    Pattern for input reads (e.g., 'data/*.fastq.gz')
    --reference FILE   Reference FASTA file
    --output DIR       Output directory (default: results)
    --profile PROFILE  Nextflow profile (default: docker)
    --cleanup-work     Remove work directory after completion
    --cleanup-all      Remove work and .nextflow directories after completion
    --help             Show this help message

Examples:
    # Setup environment only
    $0 --setup-only
    
    # Run test analysis
    $0 --test
    
    # Run with custom data
    $0 --reads 'data/*.fastq.gz' --reference reference.fasta --output my_results
    
    # Run with Singularity instead of Docker
    $0 --test --profile singularity

EOF
                exit 0
                ;;
            *)
                print_error "Unknown option: $1"
                print_error "Use --help for usage information"
                exit 1
                ;;
        esac
    done
    
    # Check system requirements
    check_requirements
    
    # Setup environment
    if ! command_exists conda; then
        install_conda
    fi
    
    create_environment
    setup_directories
    
    # Copy pipeline files to scripts directory
    if [ -f "variant_analyzer.py" ]; then
        cp variant_analyzer.py scripts/
        print_status "Copied variant analyzer to scripts/"
    fi
    
    if $SETUP_ONLY; then
        print_status "Setup completed. Environment ready for pipeline execution."
        print_status "Activate environment with: conda activate $CONDA_ENV_NAME"
        exit 0
    fi
    
    # Setup test data if requested
    if $TEST_RUN; then
        download_test_data
        READS_PATTERN="test_data/*.fastq.gz"
        REFERENCE_FILE="test_data/reference.fasta"
        OUTPUT_DIR="test_results"
        print_status "Test run configured"
    fi
    
    # Validate inputs
    if [ -z "$READS_PATTERN" ] || [ -z "$REFERENCE_FILE" ]; then
        print_error "Missing required arguments. Use --reads and --reference, or --test for test run"
        exit 1
    fi
    
    validate_inputs "$READS_PATTERN" "$REFERENCE_FILE" || exit 1
    
    # Run the pipeline
    run_pipeline "$READS_PATTERN" "$REFERENCE_FILE" "$OUTPUT_DIR" "$PROFILE" || exit 1
    
    # Run post-analysis
    run_post_analysis "$OUTPUT_DIR" "$REFERENCE_FILE"
    
    # Generate final report
    generate_final_report "$OUTPUT_DIR"
    
    # Cleanup if requested
    if $CLEANUP_WORK || $CLEANUP_NEXTFLOW; then
        cleanup
    fi
    
    print_header "=== Pipeline Execution Complete ==="
    print_status "Results available in: $OUTPUT_DIR"
    print_status "Open $OUTPUT_DIR/10_final_report/final_report.html to view results"
    
    # Display quick summary
    if [ -f "$OUTPUT_DIR/10_final_report/summary_statistics.txt" ]; then
        print_header "Quick Summary:"
        cat "$OUTPUT_DIR/10_final_report/summary_statistics.txt"
    fi
}

# Handle script interruption
trap 'print_error "Script interrupted"; exit 1' INT TERM

# Check if running as main script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi