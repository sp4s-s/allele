use anyhow::{Context, Result};
use chrono::Local;
use csv::Writer;
use serde_json::to_string_pretty;
use std::fs;
use std::path::Path;

use crate::analysis::AnalysisResults;
use crate::types::*;

/// Supported report formats
#[derive(Debug, Clone, Copy)]
pub enum ReportFormat {
    Html,
    Pdf,
    Csv,
    Json,
    Tsv,
    All,
}

/// Report generator for analysis results
pub struct ReportGenerator {
    output_dir: String,
}

impl ReportGenerator {
    pub fn new(output_dir: &Path) -> Self {
        // Create output directory if it doesn't exist
        if !output_dir.exists() {
            fs::create_dir_all(output_dir).expect("Failed to create output directory");
        }

        Self {
            output_dir: output_dir.to_string_lossy().to_string(),
        }
    }

    /// Generate reports in specified format(s)
    pub fn generate(&self, results: &AnalysisResults, format: ReportFormat) -> Result<()> {
        match format {
            ReportFormat::Html => self.generate_html_report(results)?,
            ReportFormat::Pdf => self.generate_pdf_report(results)?,
            ReportFormat::Csv => self.generate_csv_report(results)?,
            ReportFormat::Json => self.generate_json_report(results)?,
            ReportFormat::Tsv => self.generate_tsv_report(results)?,
            ReportFormat::All => {
                self.generate_html_report(results)?;
                self.generate_pdf_report(results)?;
                self.generate_csv_report(results)?;
                self.generate_json_report(results)?;
                self.generate_tsv_report(results)?;
            }
        }

        Ok(())
    }

    fn generate_html_report(&self, results: &AnalysisResults) -> Result<()> {
        let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();
        let filename = format!("{}/report_{}.html", self.output_dir, timestamp);

        let html_content = self.create_html_content(results);
        fs::write(&filename, html_content)
            .with_context(|| format!("Failed to write HTML report to {}", filename))?;

        Ok(())
    }

    fn create_html_content(&self, results: &AnalysisResults) -> String {
        let timestamp = Local::now().format("%Y-%m-%d %H:%M:%S").to_string();

        format!(
            r#"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Genetic Analysis Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 40px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 0 10px rgba(0,0,0,0.1);
        }}
        h1, h2, h3 {{
            color: #2c3e50;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }}
        th {{
            background-color: #3498db;
            color: white;
        }}
        tr:nth-child(even) {{
            background-color: #f2f2f2;
        }}
        .section {{
            margin: 30px 0;
        }}
        .summary-box {{
            background-color: #e8f4f8;
            padding: 20px;
            border-radius: 5px;
            margin: 20px 0;
        }}
        .compatibility-high {{
            background-color: #d4edda;
        }}
        .compatibility-medium {{
            background-color: #fff3cd;
        }}
        .compatibility-low {{
            background-color: #f8d7da;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Genetic Analysis Report</h1>
        <p>Generated on: {}</p>
        
        <div class="summary-box">
            <h2>Summary</h2>
            <p>This report contains genetic compatibility analysis results for {} samples.</p>
        </div>
        
        {}
        {}
        {}
        {}
        {}
    </div>
</body>
</html>"#,
            timestamp,
            results.organ_compatibility.len(),
            self.generate_organ_compatibility_html(&results.organ_compatibility),
            self.generate_disease_risk_html(&results.disease_risks),
            self.generate_hla_matching_html(&results.hla_matches),
            self.generate_ibd_detection_html(&results.ibd_segments),
            self.generate_pharmacogenomics_html(&results.pharmacogenomics)
        )
    }

    fn generate_organ_compatibility_html(&self, results: &[OrganCompatibilityResult]) -> String {
        if results.is_empty() {
            return "<div class=\"section\"><h2>Organ Compatibility</h2><p>No results available.</p></div>".to_string();
        }

        let mut html = "<div class=\"section\"><h2>Organ Compatibility</h2>\n<table>\n<tr><th>Sample</th><th>Organ</th><th>Compatibility Score</th><th>HLA Matches</th><th>Status</th></tr>\n".to_string();

        for result in results {
            let status_class = if result.compatibility_score > 0.8 {
                "compatibility-high"
            } else if result.compatibility_score > 0.6 {
                "compatibility-medium"
            } else {
                "compatibility-low"
            };

            html.push_str(&format!(
                "<tr class=\"{}\"><td>{}</td><td>{}</td><td>{:.2}%</td><td>{}/{}</td><td>{:?}</td></tr>\n",
                status_class,
                result.sample_id,
                result.organ,
                result.compatibility_score * 100.0,
                result.hla_matches.total_matches,
                result.hla_matches.total_matches + result.hla_matches.mismatch_count,
                result.crossmatch_result
            ));
        }

        html.push_str("</table>\n</div>\n");
        html
    }

    fn generate_disease_risk_html(&self, results: &[DiseaseRiskResult]) -> String {
        if results.is_empty() {
            return "<div class=\"section\"><h2>Disease Risk Analysis</h2><p>No results available.</p></div>".to_string();
        }

        let mut html = "<div class=\"section\"><h2>Disease Risk Analysis</h2>\n<table>\n<tr><th>Sample</th><th>Disease</th><th>Risk Level</th><th>Risk Score</th></tr>\n".to_string();

        for result in results {
            html.push_str(&format!(
                "<tr><td>{}</td><td>{}</td><td>{:?}</td><td>{:.2}</td></tr>\n",
                result.sample_id, result.disease, result.risk_level, result.risk_score
            ));
        }

        html.push_str("</table>\n</div>\n");
        html
    }

    fn generate_hla_matching_html(&self, results: &[HLAMatchResult]) -> String {
        if results.is_empty() {
            return "<div class=\"section\"><h2>HLA Matching</h2><p>No results available.</p></div>".to_string();
        }

        let mut html = "<div class=\"section\"><h2>HLA Matching</h2>\n<table>\n<tr><th>Sample</th><th>Match Score</th><th>HLA-A</th><th>HLA-B</th><th>HLA-DRB1</th></tr>\n".to_string();

        for result in results {
            html.push_str(&format!(
                "<tr><td>{}</td><td>{:.2}%</td><td>{}</td><td>{}</td><td>{}</td></tr>\n",
                result.sample_id,
                result.match_score * 100.0,
                result.hla_a.len(),
                result.hla_b.len(),
                result.hla_drb1.len()
            ));
        }

        html.push_str("</table>\n</div>\n");
        html
    }

    fn generate_ibd_detection_html(&self, results: &[IBDSegmentResult]) -> String {
        if results.is_empty() {
            return "<div class=\"section\"><h2>IBD Detection</h2><p>No results available.</p></div>".to_string();
        }

        let mut html = "<div class=\"section\"><h2>IBD Detection</h2>\n<table>\n<tr><th>Sample</th><th>Total Shared (cM)</th><th>Segments</th><th>Largest Segment</th><th>Predicted Relationship</th></tr>\n".to_string();

        for result in results {
            html.push_str(&format!(
                "<tr><td>{}</td><td>{:.2}</td><td>{}</td><td>{:.2} cM</td><td>{:?}</td></tr>\n",
                result.sample_id,
                result.total_shared_cm,
                result.segment_count,
                result.largest_segment_cm,
                result.predicted_relationship
            ));
        }

        html.push_str("</table>\n</div>\n");
        html
    }

    fn generate_pharmacogenomics_html(&self, results: &[PharmacogenomicResult]) -> String {
        if results.is_empty() {
            return "<div class=\"section\"><h2>Pharmacogenomics</h2><p>No results available.</p></div>".to_string();
        }

        let mut html = "<div class=\"section\"><h2>Pharmacogenomics</h2>\n<table>\n<tr><th>Sample</th><th>Gene</th><th>Phenotype</th><th>Activity Score</th><th>Drugs Affected</th></tr>\n".to_string();

        for result in results {
            let drugs: Vec<String> = result
                .clinical_annotations
                .iter()
                .map(|ca| ca.drug.clone())
                .collect();
            html.push_str(&format!(
                "<tr><td>{}</td><td>{}</td><td>{}</td><td>{:.2}</td><td>{}</td></tr>\n",
                result.sample_id,
                result.gene,
                result.phenotype,
                result.activity_score,
                drugs.join(", ")
            ));
        }

        html.push_str("</table>\n</div>\n");
        html
    }

    fn generate_pdf_report(&self, _results: &AnalysisResults) -> Result<()> {
        // PDF generation would require additional dependencies
        // For now, we'll create a placeholder
        let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();
        let filename = format!("{}/report_{}.pdf", self.output_dir, timestamp);

        let pdf_content = "PDF Report Generation Placeholder\n\nThis feature requires additional dependencies for PDF generation.";
        fs::write(&filename, pdf_content)
            .with_context(|| format!("Failed to write PDF report to {}", filename))?;

        Ok(())
    }

    fn generate_csv_report(&self, results: &AnalysisResults) -> Result<()> {
        self.generate_organ_compatibility_csv(&results.organ_compatibility)?;
        self.generate_disease_risk_csv(&results.disease_risks)?;
        self.generate_hla_matching_csv(&results.hla_matches)?;
        self.generate_ibd_detection_csv(&results.ibd_segments)?;
        self.generate_pharmacogenomics_csv(&results.pharmacogenomics)?;

        Ok(())
    }

    fn generate_organ_compatibility_csv(&self, results: &[OrganCompatibilityResult]) -> Result<()> {
        if results.is_empty() {
            return Ok(());
        }

        let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();
        let filename = format!("{}/organ_compatibility_{}.csv", self.output_dir, timestamp);

        let mut wtr = Writer::from_path(&filename)
            .with_context(|| format!("Failed to create CSV writer for {}", filename))?;

        // Write header
        wtr.write_record(&[
            "sample_id",
            "organ",
            "compatibility_score",
            "hla_matches",
            "total_hla_loci",
            "blood_type_compatible",
            "crossmatch_result",
        ])?;

        // Write data
        for result in results {
            wtr.write_record(&[
                &result.sample_id,
                &result.organ,
                &format!("{:.4}", result.compatibility_score),
                &result.hla_matches.total_matches.to_string(),
                &(result.hla_matches.total_matches + result.hla_matches.mismatch_count).to_string(),
                &result.blood_type_compatible.to_string(),
                &format!("{:?}", result.crossmatch_result),
            ])?;
        }

        wtr.flush()?;
        Ok(())
    }

    fn generate_disease_risk_csv(&self, results: &[DiseaseRiskResult]) -> Result<()> {
        if results.is_empty() {
            return Ok(());
        }

        let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();
        let filename = format!("{}/disease_risk_{}.csv", self.output_dir, timestamp);

        let mut wtr = Writer::from_path(&filename)?;

        // Write header
        wtr.write_record(&[
            "sample_id",
            "disease",
            "risk_level",
            "risk_score",
            "variants_found",
        ])?;

        // Write data
        for result in results {
            wtr.write_record(&[
                &result.sample_id,
                &result.disease,
                &format!("{:?}", result.risk_level),
                &format!("{:.4}", result.risk_score),
                &result.variants_found.len().to_string(),
            ])?;
        }

        wtr.flush()?;
        Ok(())
    }

    fn generate_hla_matching_csv(&self, results: &[HLAMatchResult]) -> Result<()> {
        if results.is_empty() {
            return Ok(());
        }

        let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();
        let filename = format!("{}/hla_matching_{}.csv", self.output_dir, timestamp);

        let mut wtr = Writer::from_path(&filename)?;

        // Write header
        wtr.write_record(&[
            "sample_id",
            "match_score",
            "hla_a_count",
            "hla_b_count",
            "hla_c_count",
            "hla_drb1_count",
            "hla_dqb1_count",
            "hla_dpb1_count",
        ])?;

        // Write data
        for result in results {
            wtr.write_record(&[
                &result.sample_id,
                &format!("{:.4}", result.match_score),
                &result.hla_a.len().to_string(),
                &result.hla_b.len().to_string(),
                &result.hla_c.len().to_string(),
                &result.hla_drb1.len().to_string(),
                &result.hla_dqb1.len().to_string(),
                &result.hla_dpb1.len().to_string(),
            ])?;
        }

        wtr.flush()?;
        Ok(())
    }

    fn generate_ibd_detection_csv(&self, results: &[IBDSegmentResult]) -> Result<()> {
        if results.is_empty() {
            return Ok(());
        }

        let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();
        let filename = format!("{}/ibd_detection_{}.csv", self.output_dir, timestamp);

        let mut wtr = Writer::from_path(&filename)?;

        // Write header
        wtr.write_record(&[
            "sample_id",
            "total_shared_cm",
            "segment_count",
            "largest_segment_cm",
            "predicted_relationship",
            "confidence",
        ])?;

        // Write data
        for result in results {
            wtr.write_record(&[
                &result.sample_id,
                &format!("{:.4}", result.total_shared_cm),
                &result.segment_count.to_string(),
                &format!("{:.4}", result.largest_segment_cm),
                &format!("{:?}", result.predicted_relationship),
                &format!("{:.4}", result.confidence),
            ])?;
        }

        wtr.flush()?;
        Ok(())
    }

    fn generate_pharmacogenomics_csv(&self, results: &[PharmacogenomicResult]) -> Result<()> {
        if results.is_empty() {
            return Ok(());
        }

        let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();
        let filename = format!("{}/pharmacogenomics_{}.csv", self.output_dir, timestamp);

        let mut wtr = Writer::from_path(&filename)?;

        // Write header
        wtr.write_record(&[
            "sample_id",
            "gene",
            "diplotype",
            "phenotype",
            "activity_score",
            "drugs_affected",
        ])?;

        // Write data
        for result in results {
            let drugs: Vec<String> = result
                .clinical_annotations
                .iter()
                .map(|ca| ca.drug.clone())
                .collect();
            wtr.write_record(&[
                &result.sample_id,
                &result.gene,
                &result.diplotype,
                &result.phenotype,
                &format!("{:.4}", result.activity_score),
                &drugs.join("; "),
            ])?;
        }

        wtr.flush()?;
        Ok(())
    }

    fn generate_json_report(&self, results: &AnalysisResults) -> Result<()> {
        let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();
        let filename = format!("{}/report_{}.json", self.output_dir, timestamp);

        let json_content =
            to_string_pretty(results).with_context(|| "Failed to serialize results to JSON")?;

        fs::write(&filename, json_content)
            .with_context(|| format!("Failed to write JSON report to {}", filename))?;

        Ok(())
    }

    fn generate_tsv_report(&self, results: &AnalysisResults) -> Result<()> {
        self.generate_organ_compatibility_tsv(&results.organ_compatibility)?;
        self.generate_disease_risk_tsv(&results.disease_risks)?;
        self.generate_hla_matching_tsv(&results.hla_matches)?;
        self.generate_ibd_detection_tsv(&results.ibd_segments)?;
        self.generate_pharmacogenomics_tsv(&results.pharmacogenomics)?;

        Ok(())
    }

    fn generate_organ_compatibility_tsv(&self, results: &[OrganCompatibilityResult]) -> Result<()> {
        if results.is_empty() {
            return Ok(());
        }

        let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();
        let filename = format!("{}/organ_compatibility_{}.tsv", self.output_dir, timestamp);

        let mut wtr = Writer::from_path(&filename)?;
        wtr.write_record(&[
            "sample_id",
            "organ",
            "compatibility_score",
            "hla_matches",
            "total_hla_loci",
            "blood_type_compatible",
            "crossmatch_result",
        ])?;

        for result in results {
            wtr.write_record(&[
                &result.sample_id,
                &result.organ,
                &format!("{:.4}", result.compatibility_score),
                &result.hla_matches.total_matches.to_string(),
                &(result.hla_matches.total_matches + result.hla_matches.mismatch_count).to_string(),
                &result.blood_type_compatible.to_string(),
                &format!("{:?}", result.crossmatch_result),
            ])?;
        }

        wtr.flush()?;
        Ok(())
    }

    fn generate_disease_risk_tsv(&self, results: &[DiseaseRiskResult]) -> Result<()> {
        if results.is_empty() {
            return Ok(());
        }

        let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();
        let filename = format!("{}/disease_risk_{}.tsv", self.output_dir, timestamp);

        let mut wtr = Writer::from_path(&filename)?;
        wtr.write_record(&[
            "sample_id",
            "disease",
            "risk_level",
            "risk_score",
            "variants_found",
        ])?;

        for result in results {
            wtr.write_record(&[
                &result.sample_id,
                &result.disease,
                &format!("{:?}", result.risk_level),
                &format!("{:.4}", result.risk_score),
                &result.variants_found.len().to_string(),
            ])?;
        }

        wtr.flush()?;
        Ok(())
    }

    fn generate_hla_matching_tsv(&self, results: &[HLAMatchResult]) -> Result<()> {
        if results.is_empty() {
            return Ok(());
        }

        let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();
        let filename = format!("{}/hla_matching_{}.tsv", self.output_dir, timestamp);

        let mut wtr = Writer::from_path(&filename)?;
        wtr.write_record(&[
            "sample_id",
            "match_score",
            "hla_a_count",
            "hla_b_count",
            "hla_c_count",
            "hla_drb1_count",
            "hla_dqb1_count",
            "hla_dpb1_count",
        ])?;

        for result in results {
            wtr.write_record(&[
                &result.sample_id,
                &format!("{:.4}", result.match_score),
                &result.hla_a.len().to_string(),
                &result.hla_b.len().to_string(),
                &result.hla_c.len().to_string(),
                &result.hla_drb1.len().to_string(),
                &result.hla_dqb1.len().to_string(),
                &result.hla_dpb1.len().to_string(),
            ])?;
        }

        wtr.flush()?;
        Ok(())
    }

    fn generate_ibd_detection_tsv(&self, results: &[IBDSegmentResult]) -> Result<()> {
        if results.is_empty() {
            return Ok(());
        }

        let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();
        let filename = format!("{}/ibd_detection_{}.tsv", self.output_dir, timestamp);

        let mut wtr = Writer::from_path(&filename)?;
        wtr.write_record(&[
            "sample_id",
            "total_shared_cm",
            "segment_count",
            "largest_segment_cm",
            "predicted_relationship",
            "confidence",
        ])?;

        for result in results {
            wtr.write_record(&[
                &result.sample_id,
                &format!("{:.4}", result.total_shared_cm),
                &result.segment_count.to_string(),
                &format!("{:.4}", result.largest_segment_cm),
                &format!("{:?}", result.predicted_relationship),
                &format!("{:.4}", result.confidence),
            ])?;
        }

        wtr.flush()?;
        Ok(())
    }

    fn generate_pharmacogenomics_tsv(&self, results: &[PharmacogenomicResult]) -> Result<()> {
        if results.is_empty() {
            return Ok(());
        }

        let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();
        let filename = format!("{}/pharmacogenomics_{}.tsv", self.output_dir, timestamp);

        let mut wtr = Writer::from_path(&filename)?;
        wtr.write_record(&[
            "sample_id",
            "gene",
            "diplotype",
            "phenotype",
            "activity_score",
            "drugs_affected",
        ])?;

        for result in results {
            let drugs: Vec<String> = result
                .clinical_annotations
                .iter()
                .map(|ca| ca.drug.clone())
                .collect();
            wtr.write_record(&[
                &result.sample_id,
                &result.gene,
                &result.diplotype,
                &result.phenotype,
                &format!("{:.4}", result.activity_score),
                &drugs.join("; "),
            ])?;
        }

        wtr.flush()?;
        Ok(())
    }
}
