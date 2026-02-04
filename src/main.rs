use anyhow::Result;
use clap::{CommandFactory, Parser, Subcommand, ValueHint};
use clap_complete::{generate, Shell};
use console::style;
use dialoguer::{theme::ColorfulTheme, Confirm, Input, Select};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::io;
use std::path::PathBuf;
use tracing::{info, warn};

mod analysis;
mod discovery;
mod output;
mod parsers;
mod types;

use analysis::{
    AlleleComparator, AnalysisResults, DiseaseAnalyzer, HLAAnalyzer, IBDDetector,
    OrganCompatibilityChecker,
};
use discovery::FileDiscovery;
use output::ReportGenerator;
use parsers::{FileParser, ParsedGeneticData};
use types::OrganType;

/// Multithreaded allele compatibility and discovery tool
#[derive(Parser, Debug)]
#[command(
    name = "allele-compatibility",
    version,
    about = "Fast genetic allele compatibility analysis with multithreading support",
    long_about = r#"
A comprehensive tool for:
- Allele compatibility checking between genetic samples
- Organ transplant compatibility analysis
- Genetic disease detection
- HLA typing and matching
- IBD (Identity by Descent) detection
- Pharmacogenomic variant analysis (PharmCAT integration)

Supports 10+ file formats including VCF, BED, FASTA, 23andMe, AncestryDNA, and more.
"#
)]
#[command(arg_required_else_help = true)]
struct Cli {
    /// Patient genetic data file
    #[arg(short, long, value_name = "FILE", value_hint = ValueHint::FilePath)]
    patient: Option<PathBuf>,

    /// Comparison files or directories (space-separated, supports tab completion)
    #[arg(short, long, value_name = "FILES", num_args = 1.., value_hint = ValueHint::AnyPath)]
    compare: Vec<PathBuf>,

    /// Recursive search for genetic data files
    #[arg(short, long, help = "Recursively search directories")]
    recursive: bool,

    /// Interactive mode with prompts for all parameters
    #[arg(short, long, help = "Interactive mode with default values")]
    interactive: bool,

    /// Number of threads (0 = auto-detect)
    #[arg(
        short,
        long,
        default_value = "0",
        help = "Number of threads (0 = auto)"
    )]
    threads: usize,

    /// Output format
    #[arg(short, long, value_enum, default_value = "html")]
    format: OutputFormat,

    /// Output directory for reports
    #[arg(short, long, default_value = "./reports")]
    output: PathBuf,

    /// Analysis type
    #[arg(short = 'a', long, value_enum, default_value = "all")]
    analysis: AnalysisType,

    /// Minimum shared DNA in centimorgans for relationship detection
    #[arg(long, default_value = "7.0")]
    min_cm: f64,

    /// Minimum segment length for IBD detection
    #[arg(long, default_value = "500")]
    min_segment_length: usize,

    /// Organ type for transplant compatibility
    #[arg(long, value_enum)]
    organ: Option<OrganType>,

    /// Enable verbose logging
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,

    /// Generate shell completions
    #[arg(long, value_enum, value_name = "SHELL")]
    completions: Option<Shell>,

    /// Subcommands
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Generate shell completions
    Completions { shell: Shell },
    /// List supported file formats
    Formats,
    /// Check external tool dependencies
    CheckDeps,
}

#[derive(Debug, Clone, Copy, Default, clap::ValueEnum)]
enum OutputFormat {
    #[default]
    Html,
    Pdf,
    Csv,
    Json,
    Tsv,
    All,
}

impl From<OutputFormat> for output::ReportFormat {
    fn from(format: OutputFormat) -> output::ReportFormat {
        match format {
            OutputFormat::Html => output::ReportFormat::Html,
            OutputFormat::Pdf => output::ReportFormat::Pdf,
            OutputFormat::Csv => output::ReportFormat::Csv,
            OutputFormat::Json => output::ReportFormat::Json,
            OutputFormat::Tsv => output::ReportFormat::Tsv,
            OutputFormat::All => output::ReportFormat::All,
        }
    }
}

#[derive(Debug, Clone, Copy, Default, clap::ValueEnum)]
enum AnalysisType {
    #[default]
    All,
    OrganCompatibility,
    Disease,
    HLA,
    IBD,
    Pharmacogenomics,
    Relationship,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Handle shell completions
    if let Some(shell) = cli.completions {
        generate_completions(shell);
        return Ok(());
    }

    if let Some(Commands::Completions { shell }) = cli.command {
        generate_completions(shell);
        return Ok(());
    }

    if let Some(Commands::Formats) = cli.command {
        list_formats();
        return Ok(());
    }

    if let Some(Commands::CheckDeps) = cli.command {
        check_dependencies()?;
        return Ok(());
    }

    // Initialize logging
    init_logging(cli.verbose);

    // Run interactive mode if requested
    let config = if cli.interactive {
        run_interactive_mode()?
    } else {
        AppConfig::from_cli(&cli)
    };

    // Initialize thread pool
    init_thread_pool(config.threads)?;

    info!("Starting allele compatibility analysis...");
    info!("Using {} threads", rayon::current_num_threads());

    // Run the analysis
    run_analysis(config)?;

    Ok(())
}

fn generate_completions(shell: Shell) {
    let mut cmd = Cli::command();
    let name = cmd.get_name().to_string();
    generate(shell, &mut cmd, name, &mut io::stdout());
}

fn list_formats() {
    println!("{}", style("Supported Genetic Data Formats:").bold().cyan());
    println!();

    let formats = vec![
        (
            "VCF",
            "Variant Call Format (.vcf, .vcf.gz)",
            "Standard genomic variant format",
        ),
        ("BCF", "Binary Call Format (.bcf)", "Binary compressed VCF"),
        (
            "BED",
            "Browser Extensible Data (.bed)",
            "Genomic intervals and annotations",
        ),
        (
            "FASTA",
            "FASTA (.fa, .fasta, .fna)",
            "Raw DNA/RNA sequences",
        ),
        (
            "FASTQ",
            "FASTQ (.fq, .fastq)",
            "Sequencing reads with quality scores",
        ),
        (
            "23andMe",
            "23andMe Raw Data (.txt)",
            "Direct-to-consumer genetic data",
        ),
        (
            "AncestryDNA",
            "AncestryDNA (.txt, .csv)",
            "Ancestry.com genetic data",
        ),
        (
            "FTDNA",
            "Family Tree DNA (.csv)",
            "FTDNA autosomal transfer",
        ),
        ("MyHeritage", "MyHeritage (.csv)", "MyHeritage DNA raw data"),
        (
            "PLINK",
            "PLINK (.ped, .map, .bed, .bim, .fam)",
            "GWAS format",
        ),
        (
            "GVF",
            "Genome Variation Format (.gvf)",
            "Variant annotation format",
        ),
        (
            "GFF",
            "General Feature Format (.gff, .gff3)",
            "Gene annotations",
        ),
        (
            "SAM/BAM",
            "Sequence Alignment (.sam, .bam, .cram)",
            "Read alignments",
        ),
        ("GenBank", "GenBank (.gb, .gbk)", "Annotated sequences"),
    ];

    for (name, ext, desc) in formats {
        println!("  {} - {}", style(name).green().bold(), style(ext).yellow());
        println!("         {}", style(desc).dim());
    }
}

fn check_dependencies() -> Result<()> {
    println!("{}", style("Checking External Dependencies:").bold().cyan());
    println!();

    let tools = vec![
        (
            "PharmCAT",
            "pharmcat",
            vec!["pharmcat", "java -jar pharmcat"],
        ),
        ("Athlates", "athlates", vec!["athlates"]),
        ("Seq2HLA", "seq2HLA", vec!["seq2HLA"]),
        ("PLINK", "plink", vec!["plink", "plink2"]),
        ("BCFtools", "bcftools", vec!["bcftools"]),
        ("SAMtools", "samtools", vec!["samtools"]),
        ("BEDtools", "bedtools", vec!["bedtools"]),
    ];

    for (name, cmd, alternatives) in tools {
        let found = alternatives.iter().any(|alt| {
            std::process::Command::new("sh")
                .arg("-c")
                .arg(format!(
                    "which {} 2>/dev/null || command -v {} 2>/dev/null",
                    alt, alt
                ))
                .output()
                .map(|o| o.status.success())
                .unwrap_or(false)
        });

        if found {
            println!("  {} {}", style("✓").green(), style(name).green());
        } else {
            println!("  {} {} (optional)", style("✗").red(), style(name).dim());
        }
    }

    Ok(())
}

fn init_logging(verbose: u8) {
    let level = match verbose {
        0 => "warn",
        1 => "info",
        2 => "debug",
        _ => "trace",
    };

    tracing_subscriber::fmt()
        .with_env_filter(format!("allele_compatibility={}", level))
        .init();
}

fn init_thread_pool(threads: usize) -> Result<()> {
    let num_threads = if threads == 0 {
        num_cpus::get()
    } else {
        threads
    };

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .map_err(|e| anyhow::anyhow!("Failed to initialize thread pool: {}", e))?;

    Ok(())
}

fn run_interactive_mode() -> Result<AppConfig> {
    println!(
        "{}",
        style("╔══════════════════════════════════════════════════════════════╗").cyan()
    );
    println!(
        "{}",
        style("║     Allele Compatibility Analysis - Interactive Mode         ║")
            .cyan()
            .bold()
    );
    println!(
        "{}",
        style("╚══════════════════════════════════════════════════════════════╝").cyan()
    );
    println!();

    let theme = ColorfulTheme::default();

    // Patient file selection
    let patient: String = Input::with_theme(&theme)
        .with_prompt("Patient genetic data file")
        .interact_text()?;

    // Comparison files
    let compare_input: String = Input::with_theme(&theme)
        .with_prompt("Comparison files/directories (space-separated, supports tab completion)")
        .allow_empty(true)
        .interact_text()?;

    let compare: Vec<PathBuf> = if compare_input.is_empty() {
        vec![]
    } else {
        compare_input
            .split_whitespace()
            .map(PathBuf::from)
            .collect()
    };

    // Recursive search
    let recursive = Confirm::with_theme(&theme)
        .with_prompt("Enable recursive directory search?")
        .default(true)
        .interact()?;

    // Analysis type
    let analysis_types = vec![
        "All analyses",
        "Organ transplant compatibility",
        "Genetic disease detection",
        "HLA typing",
        "IBD detection",
        "Pharmacogenomics",
        "Relationship analysis",
    ];

    let analysis_idx = Select::with_theme(&theme)
        .with_prompt("Select analysis type")
        .default(0)
        .items(&analysis_types)
        .interact()?;

    let analysis = match analysis_idx {
        0 => AnalysisType::All,
        1 => AnalysisType::OrganCompatibility,
        2 => AnalysisType::Disease,
        3 => AnalysisType::HLA,
        4 => AnalysisType::IBD,
        5 => AnalysisType::Pharmacogenomics,
        6 => AnalysisType::Relationship,
        _ => AnalysisType::All,
    };

    // Organ type if applicable
    let organ = if matches!(
        analysis,
        AnalysisType::OrganCompatibility | AnalysisType::All
    ) {
        let organs = vec![
            "Any/All organs",
            "Kidney",
            "Liver",
            "Heart",
            "Lung",
            "Pancreas",
            "Bone Marrow",
            "Cornea",
            "Skin",
        ];

        let organ_idx = Select::with_theme(&theme)
            .with_prompt("Select organ type for compatibility")
            .default(0)
            .items(&organs)
            .interact()?;

        match organ_idx {
            0 => None,
            1 => Some(OrganType::Kidney),
            2 => Some(OrganType::Liver),
            3 => Some(OrganType::Heart),
            4 => Some(OrganType::Lung),
            5 => Some(OrganType::Pancreas),
            6 => Some(OrganType::BoneMarrow),
            7 => Some(OrganType::Cornea),
            8 => Some(OrganType::Skin),
            _ => None,
        }
    } else {
        None
    };

    // Output format
    let formats = vec!["HTML", "PDF", "CSV", "JSON", "TSV", "All formats"];
    let format_idx = Select::with_theme(&theme)
        .with_prompt("Select output format")
        .default(0)
        .items(&formats)
        .interact()?;

    let format = match format_idx {
        0 => OutputFormat::Html,
        1 => OutputFormat::Pdf,
        2 => OutputFormat::Csv,
        3 => OutputFormat::Json,
        4 => OutputFormat::Tsv,
        5 => OutputFormat::All,
        _ => OutputFormat::Html,
    };

    // Output directory
    let output: String = Input::with_theme(&theme)
        .with_prompt("Output directory")
        .default("./reports".to_string())
        .interact_text()?;

    // Threads
    let threads: usize = Input::with_theme(&theme)
        .with_prompt("Number of threads (0 = auto-detect)")
        .default(0)
        .interact_text()?;

    Ok(AppConfig {
        patient: PathBuf::from(patient),
        compare,
        recursive,
        threads,
        format,
        output: PathBuf::from(output),
        analysis,
        organ,
        min_cm: 7.0,
        min_segment_length: 500,
    })
}

fn run_analysis(config: AppConfig) -> Result<()> {
    let pb = ProgressBar::new(100);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}% {msg}")?
            .progress_chars("#>-"),
    );

    // Step 1: Discover files
    pb.set_message("Discovering genetic data files...");
    let discovery = FileDiscovery::new(config.recursive);
    let files_to_process = discovery.discover(&config.patient, &config.compare)?;
    pb.set_position(10);

    info!("Found {} files to analyze", files_to_process.len());

    // Step 2: Parse patient data
    pb.set_message("Parsing patient genetic data...");
    let parser = FileParser::new();
    let patient_data = parser.parse(&config.patient)?;
    pb.set_position(20);

    // Step 3: Parse comparison files in parallel
    pb.set_message("Parsing comparison files...");
    let comparison_data: Vec<ParsedGeneticData> = files_to_process
        .par_iter()
        .filter_map(|path| match parser.parse(path) {
            Ok(data) => Some(data),
            Err(e) => {
                warn!("Failed to parse {}: {}", path.display(), e);
                None
            }
        })
        .collect();
    pb.set_position(40);

    info!(
        "Successfully parsed {} comparison files",
        comparison_data.len()
    );

    // Step 4: Run analyses based on selected type
    let mut results = AnalysisResults::new();

    match config.analysis {
        AnalysisType::All | AnalysisType::OrganCompatibility => {
            pb.set_message("Analyzing organ transplant compatibility...");
            let checker = OrganCompatibilityChecker::new(config.organ);
            results.organ_compatibility = checker.check(&patient_data, &comparison_data)?;
            pb.set_position(50);
        }
        _ => {}
    }

    match config.analysis {
        AnalysisType::All | AnalysisType::Disease => {
            pb.set_message("Analyzing genetic disease risks...");
            let analyzer = DiseaseAnalyzer::new();
            results.disease_risks = analyzer.analyze(&patient_data, &comparison_data)?;
            pb.set_position(60);
        }
        _ => {}
    }

    match config.analysis {
        AnalysisType::All | AnalysisType::HLA => {
            pb.set_message("Performing HLA typing analysis...");
            let hla = HLAAnalyzer::new();
            results.hla_matches = hla.analyze(&patient_data, &comparison_data)?;
            pb.set_position(70);
        }
        _ => {}
    }

    match config.analysis {
        AnalysisType::All | AnalysisType::IBD | AnalysisType::Relationship => {
            pb.set_message("Detecting IBD segments...");
            let ibd = IBDDetector::new(config.min_cm, config.min_segment_length);
            results.ibd_segments = ibd.detect(&patient_data, &comparison_data)?;
            pb.set_position(80);
        }
        _ => {}
    }

    match config.analysis {
        AnalysisType::All | AnalysisType::Pharmacogenomics => {
            pb.set_message("Analyzing pharmacogenomic variants...");
            let comparator = AlleleComparator::new();
            results.pharmacogenomics = comparator.compare(&patient_data, &comparison_data)?;
            pb.set_position(90);
        }
        _ => {}
    }

    // Step 5: Generate reports
    pb.set_message("Generating reports...");
    let generator = ReportGenerator::new(&config.output);
    generator.generate(&results, config.format.into())?;
    pb.set_position(100);

    pb.finish_with_message("Analysis complete!");

    println!(
        "\n{} Reports saved to: {}",
        style("✓").green().bold(),
        style(config.output.display()).cyan()
    );

    Ok(())
}

#[derive(Debug)]
struct AppConfig {
    patient: PathBuf,
    compare: Vec<PathBuf>,
    recursive: bool,
    threads: usize,
    format: OutputFormat,
    output: PathBuf,
    analysis: AnalysisType,
    organ: Option<OrganType>,
    min_cm: f64,
    min_segment_length: usize,
}

impl AppConfig {
    fn from_cli(cli: &Cli) -> Self {
        Self {
            patient: cli.patient.clone().unwrap_or_default(),
            compare: cli.compare.clone(),
            recursive: cli.recursive,
            threads: cli.threads,
            format: cli.format,
            output: cli.output.clone(),
            analysis: cli.analysis,
            organ: cli.organ,
            min_cm: cli.min_cm,
            min_segment_length: cli.min_segment_length,
        }
    }
}
