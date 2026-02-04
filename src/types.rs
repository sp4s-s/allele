use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Organ types for transplant compatibility
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, clap::ValueEnum)]
pub enum OrganType {
    Kidney,
    Liver,
    Heart,
    Lung,
    Pancreas,
    BoneMarrow,
    Cornea,
    Skin,
}

/// Represents a single genetic variant/SNP
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Variant {
    pub chromosome: String,
    pub position: u64,
    pub rsid: Option<String>,
    pub reference_allele: String,
    pub alternate_alleles: Vec<String>,
    pub genotype: Genotype,
    pub quality: Option<f32>,
    pub filter: Option<String>,
    pub info: HashMap<String, String>,
}

/// Genotype representation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Genotype {
    HomozygousReference, // AA (ref/ref)
    HomozygousAlternate, // BB (alt/alt)
    Heterozygous,        // AB (ref/alt)
    HeterozygousAltAlt,  // AC (alt1/alt2) - multi-allelic
    NoCall,              // ./. or missing
    Partial,             // . or incomplete
}

impl Genotype {
    pub fn from_string(s: &str) -> Self {
        match s {
            "0/0" | "0|0" | "AA" => Genotype::HomozygousReference,
            "1/1" | "1|1" | "BB" => Genotype::HomozygousAlternate,
            "0/1" | "0|1" | "1/0" | "1|0" | "AB" | "BA" => Genotype::Heterozygous,
            "./." | ".|." | "" => Genotype::NoCall,
            _ => {
                // Handle other cases like 0/2, 1/2, etc.
                if s.contains('/') || s.contains('|') {
                    let parts: Vec<&str> = s.split(&['/', '|'][..]).collect();
                    if parts.len() == 2 {
                        if parts[0] == parts[1] {
                            if parts[0] == "0" {
                                Genotype::HomozygousReference
                            } else {
                                Genotype::HomozygousAlternate
                            }
                        } else if parts[0] == "0" || parts[1] == "0" {
                            Genotype::Heterozygous
                        } else {
                            Genotype::HeterozygousAltAlt
                        }
                    } else {
                        Genotype::NoCall
                    }
                } else {
                    Genotype::Partial
                }
            }
        }
    }

    pub fn is_homozygous(&self) -> bool {
        matches!(
            self,
            Genotype::HomozygousReference | Genotype::HomozygousAlternate
        )
    }

    pub fn is_heterozygous(&self) -> bool {
        matches!(self, Genotype::Heterozygous | Genotype::HeterozygousAltAlt)
    }

    pub fn is_no_call(&self) -> bool {
        matches!(self, Genotype::NoCall | Genotype::Partial)
    }

    /// Check if two genotypes are compatible for organ transplant
    /// For most organs, we need at least one shared allele
    pub fn is_compatible(&self, other: &Genotype) -> bool {
        match (self, other) {
            // No-calls are treated as wildcards (could match anything)
            (Genotype::NoCall, _) | (_, Genotype::NoCall) => true,

            // Homozygous reference matches with anything that has reference allele
            (Genotype::HomozygousReference, Genotype::HomozygousReference) => true,
            (Genotype::HomozygousReference, Genotype::Heterozygous) => true,
            (Genotype::Heterozygous, Genotype::HomozygousReference) => true,

            // Homozygous alternate matches with anything that has alternate allele
            (Genotype::HomozygousAlternate, Genotype::HomozygousAlternate) => true,
            (Genotype::HomozygousAlternate, Genotype::Heterozygous) => true,
            (Genotype::Heterozygous, Genotype::HomozygousAlternate) => true,

            // Heterozygous matches with heterozygous (at least one shared)
            (Genotype::Heterozygous, Genotype::Heterozygous) => true,

            // Complex cases - assume compatible for safety
            _ => true,
        }
    }
}

/// HLA allele representation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HLAAllele {
    pub gene: String,   // e.g., "HLA-A", "HLA-B", "HLA-DRB1"
    pub allele: String, // e.g., "*01:01"
    pub resolution: HLAResolution,
    pub confidence: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum HLAResolution {
    TwoDigit,   // A*01
    FourDigit,  // A*01:01
    SixDigit,   // A*01:01:01
    EightDigit, // A*01:01:01:01
    Unknown,
}

/// IBD (Identity by Descent) segment
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IBDSegment {
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub length_cm: f64,
    pub snp_count: usize,
    pub start_position_bp: u64,
    pub end_position_bp: u64,
}

/// Pharmacogenomic variant (for PharmCAT integration)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PharmacoVariant {
    pub gene: String, // e.g., "CYP2C19"
    pub rsid: String,
    pub position: u64,
    pub reference: String,
    pub alternate: String,
    pub phenotype: Option<String>,
    pub activity_score: Option<f32>,
}

/// Disease risk variant
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiseaseVariant {
    pub disease: String,
    pub rsid: String,
    pub chromosome: String,
    pub position: u64,
    pub risk_allele: String,
    pub odds_ratio: Option<f32>,
    pub risk_category: RiskCategory,
    pub confidence: ConfidenceLevel,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RiskCategory {
    High,
    Moderate,
    Low,
    Protective,
    Unknown,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ConfidenceLevel {
    Definitive,
    Likely,
    Uncertain,
    Unknown,
}

/// Organ compatibility result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrganCompatibilityResult {
    pub sample_id: String,
    pub organ: String,
    pub compatibility_score: f32,
    pub hla_matches: HLATypingResult,
    pub blood_type_compatible: bool,
    pub crossmatch_result: CrossmatchResult,
    pub recommendations: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HLATypingResult {
    pub a_matches: usize,
    pub b_matches: usize,
    pub dr_matches: usize,
    pub dq_matches: usize,
    pub dp_matches: usize,
    pub total_matches: usize,
    pub mismatch_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CrossmatchResult {
    Compatible,
    Incompatible,
    RequiresFurtherTesting,
    Unknown,
}

/// Disease risk analysis result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiseaseRiskResult {
    pub sample_id: String,
    pub disease: String,
    pub risk_level: RiskLevel,
    pub risk_score: f32,
    pub variants_found: Vec<DiseaseVariant>,
    pub recommendations: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RiskLevel {
    VeryHigh,
    High,
    Moderate,
    Low,
    VeryLow,
    Unknown,
}

/// HLA match result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HLAMatchResult {
    pub sample_id: String,
    pub hla_a: Vec<HLAAllele>,
    pub hla_b: Vec<HLAAllele>,
    pub hla_c: Vec<HLAAllele>,
    pub hla_drb1: Vec<HLAAllele>,
    pub hla_dqb1: Vec<HLAAllele>,
    pub hla_dpb1: Vec<HLAAllele>,
    pub match_score: f32,
}

/// IBD detection result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IBDSegmentResult {
    pub sample_id: String,
    pub total_shared_cm: f64,
    pub segment_count: usize,
    pub largest_segment_cm: f64,
    pub segments: Vec<IBDSegment>,
    pub predicted_relationship: RelationshipPrediction,
    pub confidence: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RelationshipPrediction {
    IdenticalTwin,
    ParentChild,
    FullSibling,
    HalfSibling,
    GrandparentGrandchild,
    AuntUncleNieceNephew,
    FirstCousin,
    FirstCousinOnceRemoved,
    SecondCousin,
    SecondCousinOnceRemoved,
    ThirdCousin,
    Distant,
    Unrelated,
    Unknown,
}

impl RelationshipPrediction {
    pub fn from_shared_cm(total_cm: f64) -> Self {
        match total_cm {
            cm if cm > 3400.0 => RelationshipPrediction::IdenticalTwin,
            cm if cm > 3100.0 => RelationshipPrediction::ParentChild,
            cm if cm > 2200.0 => RelationshipPrediction::FullSibling,
            cm if cm > 1300.0 => RelationshipPrediction::HalfSibling,
            cm if cm > 1150.0 => RelationshipPrediction::GrandparentGrandchild,
            cm if cm > 800.0 => RelationshipPrediction::AuntUncleNieceNephew,
            cm if cm > 400.0 => RelationshipPrediction::FirstCousin,
            cm if cm > 200.0 => RelationshipPrediction::FirstCousinOnceRemoved,
            cm if cm > 100.0 => RelationshipPrediction::SecondCousin,
            cm if cm > 50.0 => RelationshipPrediction::SecondCousinOnceRemoved,
            cm if cm > 30.0 => RelationshipPrediction::ThirdCousin,
            cm if cm > 10.0 => RelationshipPrediction::Distant,
            _ => RelationshipPrediction::Unrelated,
        }
    }
}

/// Pharmacogenomic analysis result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PharmacogenomicResult {
    pub sample_id: String,
    pub gene: String,
    pub diplotype: String,
    pub phenotype: String,
    pub activity_score: f32,
    pub clinical_annotations: Vec<ClinicalAnnotation>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClinicalAnnotation {
    pub drug: String,
    pub implication: String,
    pub recommendation: String,
    pub evidence_level: String,
}

/// Sample metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleMetadata {
    pub sample_id: String,
    pub source_file: String,
    pub file_format: FileFormat,
    pub genome_build: String,
    pub sequencing_platform: Option<String>,
    pub call_rate: Option<f32>,
    pub heterozygosity: Option<f32>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum FileFormat {
    VCF,
    BCF,
    BED,
    FASTA,
    FASTQ,
    AndMe,
    AncestryDNA,
    FTDNA,
    MyHeritage,
    PLINK,
    GVF,
    GFF,
    SAM,
    BAM,
    CRAM,
    GenBank,
    CSV,
    TSV,
    TXT,
    Unknown,
}

impl FileFormat {
    pub fn from_extension(ext: &str) -> Self {
        match ext.to_lowercase().as_str() {
            "vcf" | "vcf.gz" => FileFormat::VCF,
            "bcf" => FileFormat::BCF,
            "bed" => FileFormat::BED,
            "fa" | "fasta" | "fna" => FileFormat::FASTA,
            "fq" | "fastq" => FileFormat::FASTQ,
            "txt" => FileFormat::TXT,
            "csv" => FileFormat::CSV,
            "tsv" => FileFormat::TSV,
            "ped" | "map" | "bim" | "fam" => FileFormat::PLINK,
            "gvf" => FileFormat::GVF,
            "gff" | "gff3" => FileFormat::GFF,
            "sam" => FileFormat::SAM,
            "bam" => FileFormat::BAM,
            "cram" => FileFormat::CRAM,
            "gb" | "gbk" => FileFormat::GenBank,
            _ => FileFormat::Unknown,
        }
    }

    pub fn is_genetic_format(&self) -> bool {
        !matches!(self, FileFormat::Unknown)
    }
}

/// Genetic coordinate
#[derive(Debug, Clone, Serialize, Deserialize, Hash, Eq, PartialEq)]
pub struct Coordinate {
    pub chromosome: String,
    pub position: u64,
}

impl Coordinate {
    pub fn new(chromosome: impl Into<String>, position: u64) -> Self {
        Self {
            chromosome: chromosome.into(),
            position,
        }
    }
}

/// Phased haplotype block
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HaplotypeBlock {
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub variants: Vec<PhasedVariant>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhasedVariant {
    pub position: u64,
    pub allele1: String,
    pub allele2: String,
    pub phase_set: String,
}

/// Discordant SNP result (for parent-child/trio analysis)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiscordantSNP {
    pub chromosome: String,
    pub position: u64,
    pub rsid: Option<String>,
    pub parent1_genotype: Genotype,
    pub parent2_genotype: Genotype,
    pub child_genotype: Genotype,
    pub is_mendelian_error: bool,
}

/// Shared segment for relationship analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SharedSegment {
    pub chromosome: String,
    pub start_bp: u64,
    pub end_bp: u64,
    pub length_cm: f64,
    pub snp_density: f64,
    pub start_snp: String,
    pub end_snp: String,
}

/// Population ancestry information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AncestryEstimate {
    pub population: String,
    pub percentage: f32,
    pub confidence_interval: (f32, f32),
}

/// Quality metrics for a sample
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityMetrics {
    pub snp_count: usize,
    pub no_call_rate: f32,
    pub heterozygosity_rate: f32,
    pub ti_tv_ratio: f32,
    pub read_depth_mean: Option<f32>,
    pub read_depth_std: Option<f32>,
}
