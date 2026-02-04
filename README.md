# Allele Compatibility Toolkit

A high-performance, multithreaded genetic analysis tool for allele compatibility checking, organ transplant compatibility analysis, and genetic disease risk assessment.

## Features

- **Ultra-fast Processing**: Utilizes all CPU cores for maximum performance
- **Multi-format Support**: Parses 10+ genetic data formats including VCF, 23andMe, AncestryDNA, FASTA, and more
- **Comprehensive Analysis**:
  - Organ transplant compatibility checking
  - Genetic disease risk detection
  - HLA typing and matching
  - IBD (Identity by Descent) detection
  - Pharmacogenomic variant analysis
  - Relationship prediction
- **Multiple Output Formats**: HTML, PDF, CSV, JSON, TSV reports
- **Interactive Mode**: User-friendly prompts with sensible defaults
- **Automatic File Discovery**: Recursively searches directories for genetic data files
- **Shell Completions**: Tab completion for file paths and commands

## Installation

### Prerequisites

- Rust 1.60 or later
- Cargo package manager

### Building from Source

```bash
git clone https://github.com/yourusername/allele-compatibility.git
cd allele-compatibility
cargo build --release
```

The compiled binary will be located at `target/release/allele-compatibility`.

## Usage

### Basic Usage

```bash
# Compare patient file with all files in a directory
./allele-compatibility -p patient.vcf -c comparison_data/

# Recursive search with HTML report output
./allele-compatibility -p patient.vcf -c data/ -r -f html -o reports/

# Interactive mode
./allele-compatibility -i
```

### Command Line Options

```
Usage: allele-compatibility [OPTIONS]

Options:
  -p, --patient <FILE>            Patient genetic data file
  -c, --compare <FILES>...        Comparison files or directories (space-separated)
  -r, --recursive                 Recursively search directories
  -i, --interactive               Interactive mode with default values
  -t, --threads <THREADS>         Number of threads (0 = auto) [default: 0]
  -f, --format <FORMAT>           Output format [default: html] [possible values: html, pdf, csv, json, tsv, all]
  -o, --output <OUTPUT>           Output directory for reports [default: ./reports]
  -a, --analysis <ANALYSIS>       Analysis type [default: all] [possible values: all, organ-compatibility, disease, hla, ibd, pharmacogenomics, relationship]
      --min-cm <MIN_CM>           Minimum shared DNA in centimorgans for relationship detection [default: 7.0]
      --min-segment-length <MIN_SEGMENT_LENGTH>  Minimum segment length for IBD detection [default: 500]
      --organ <ORGAN>             Organ type for transplant compatibility [possible values: kidney, liver, heart, lung, pancreas, bone-marrow, cornea, skin]
  -v, --verbose                   Enable verbose logging
      --completions <SHELL>       Generate shell completions [possible values: bash, elvish, fish, powershell, zsh]
  -h, --help                      Print help
  -V, --version                   Print version
```

### Subcommands

```
# List supported file formats
./allele-compatibility formats

# Check external tool dependencies
./allele-compatibility check-deps

# Generate shell completions
./allele-compatibility --completions bash > allele-compatibility.bash
```

## Supported File Formats

- **VCF/BCF**: Variant Call Format (standard genomic variant format)
- **23andMe**: Direct-to-consumer genetic data
- **AncestryDNA**: Ancestry.com genetic data
- **FTDNA**: Family Tree DNA autosomal transfer
- **MyHeritage**: MyHeritage DNA raw data
- **PLINK**: GWAS format (PED/MAP, BED/BIM/FAM)
- **FASTA/FASTQ**: Raw DNA/RNA sequences and reads
- **BED**: Browser Extensible Data (genomic intervals)
- **GVF/GFF**: Genome Variation Format and Gene annotations
- **SAM/BAM/CRAM**: Sequence alignment formats

## Analysis Modules

### Organ Transplant Compatibility

Analyzes HLA matching, blood type compatibility, and generates compatibility scores for organ transplants.

### Genetic Disease Risk

Identifies disease-associated variants and calculates risk scores for various genetic conditions.

### HLA Typing

Performs HLA allele matching for immune system compatibility analysis.

### IBD Detection

Detects Identity by Descent segments to predict relationships and shared ancestry.

### Pharmacogenomics

Analyzes drug metabolism variants to predict medication responses and optimal dosages.

## Performance

The toolkit automatically detects and utilizes all available CPU cores for parallel processing. For large datasets:
- VCF files with millions of variants are processed efficiently
- Multi-file comparisons are distributed across threads
- Memory usage is optimized through streaming parsers

## External Tool Integration

The toolkit can integrate with specialized bioinformatics tools:
- **PharmCAT**: Pharmacogenomic variant matching
- **Athlates/Seq2HLA**: HLA typing from sequencing data
- **PLINK**: Population genetics analysis

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Built with Rust for performance and memory safety
- Uses Rayon for efficient parallel processing
- Inspired by tools like PLINK, HLA*LA, and PharmCAT