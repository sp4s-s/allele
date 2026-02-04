//! # Allele Compatibility Toolkit
//!
//! A comprehensive toolkit for genetic allele compatibility analysis, organ transplant
//! compatibility checking, and genetic disease risk assessment.
//!
//! ## Features
//!
//! - Multi-threaded processing for fast analysis of large genetic datasets
//! - Support for 10+ genetic data formats (VCF, 23andMe, AncestryDNA, etc.)
//! - Organ transplant compatibility analysis
//! - Genetic disease risk detection
//! - HLA typing and matching
//! - IBD (Identity by Descent) detection
//! - Pharmacogenomic variant analysis
//! - Multiple output formats (HTML, PDF, CSV, JSON, TSV)

pub mod analysis;
pub mod discovery;
pub mod output;
pub mod parsers;
pub mod types;

// Re-export key types
pub use analysis::{
    AlleleComparator, DiseaseAnalyzer, HLAAnalyzer, IBDDetector, OrganCompatibilityChecker,
};
pub use discovery::FileDiscovery;
pub use output::{ReportFormat, ReportGenerator};
pub use parsers::{FileParser, ParsedGeneticData};
pub use types::*;