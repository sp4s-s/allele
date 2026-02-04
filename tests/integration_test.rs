use allele_compatibility::{
    analysis::*, discovery::FileDiscovery, output::*, parsers::*, types::*,
};
use std::path::PathBuf;

#[test]
fn test_project_compiles() {
    // This test just verifies that the project compiles and exports the expected types
    let _discovery = FileDiscovery::new(true);
    let _parser = FileParser::new();
    let _comparator = AlleleComparator::new();
    let _disease_analyzer = DiseaseAnalyzer::new();
    let _hla_analyzer = HLAAnalyzer::new();
    let _ibd_detector = IBDDetector::new(7.0, 500);
    let _organ_checker = OrganCompatibilityChecker::new(None);

    // Verify enums are accessible
    let _format = ReportFormat::Html;
    let _genotype = Genotype::HomozygousReference;
    let _risk_level = RiskLevel::Moderate;

    assert!(true); // If we get here, compilation succeeded
}
