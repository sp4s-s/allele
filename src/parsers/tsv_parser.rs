use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;
use std::path::Path;

use crate::parsers::{
    detect_delimiter, normalize_chromosome, open_file, parse_genotype, ParsedGeneticData,
};
use crate::types::*;

/// Generic TSV/CSV parser for genetic data
pub struct TsvParser;

impl TsvParser {
    pub fn new() -> Self {
        Self
    }

    pub fn parse(&self, path: &Path) -> Result<ParsedGeneticData> {
        let mut reader = open_file(path)?;
        let mut data = self.initialize_data(path)?;

        // Read first line to determine delimiter and column positions
        let mut first_line = String::new();
        reader.read_line(&mut first_line)?;

        let delimiter = detect_delimiter(&first_line);
        let headers: Vec<&str> = first_line.trim().split(delimiter).collect();
        let column_mapping = self.map_columns(&headers)?;

        let mut line = String::new();
        while reader.read_line(&mut line)? > 0 {
            self.parse_data_line(&line, delimiter, &column_mapping, &mut data)?;
            line.clear();
        }

        data.calculate_quality_metrics();
        Ok(data)
    }

    fn initialize_data(&self, path: &Path) -> Result<ParsedGeneticData> {
        let filename = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();

        let format = if path.extension().map(|e| e == "csv").unwrap_or(false) {
            FileFormat::CSV
        } else {
            FileFormat::TSV
        };

        Ok(ParsedGeneticData::new(
            filename.clone(),
            path.to_string_lossy().to_string(),
            format,
        ))
    }

    fn map_columns(&self, headers: &[&str]) -> Result<HashMap<String, usize>> {
        let mut mapping = HashMap::new();

        for (i, header) in headers.iter().enumerate() {
            let header_lower = header.trim().to_lowercase();

            // Map common column names
            match header_lower.as_str() {
                "chromosome" | "chr" | "chrom" => {
                    mapping.insert("chromosome".to_string(), i);
                }
                "position" | "pos" | "bp" => {
                    mapping.insert("position".to_string(), i);
                }
                "rsid" | "rs#" | "snp" => {
                    mapping.insert("rsid".to_string(), i);
                }
                "genotype" | "gt" | "allele1" | "allele2" => {
                    mapping.insert("genotype".to_string(), i);
                }
                "reference" | "ref" => {
                    mapping.insert("reference".to_string(), i);
                }
                "alternate" | "alt" => {
                    mapping.insert("alternate".to_string(), i);
                }
                _ => {}
            }
        }

        // Check that we have the minimum required columns
        if !mapping.contains_key("chromosome") || !mapping.contains_key("position") {
            return Err(anyhow!("Required columns (chromosome, position) not found"));
        }

        Ok(mapping)
    }

    fn parse_data_line(
        &self,
        line: &str,
        delimiter: char,
        column_mapping: &HashMap<String, usize>,
        data: &mut ParsedGeneticData,
    ) -> Result<()> {
        let line = line.trim();
        if line.is_empty() {
            return Ok(());
        }

        let parts: Vec<&str> = line.split(delimiter).collect();

        // Extract required fields
        let chromosome_idx = *column_mapping
            .get("chromosome")
            .ok_or_else(|| anyhow!("Chromosome column not found"))?;
        let position_idx = *column_mapping
            .get("position")
            .ok_or_else(|| anyhow!("Position column not found"))?;

        if chromosome_idx >= parts.len() || position_idx >= parts.len() {
            return Err(anyhow!("Line has insufficient columns"));
        }

        let chromosome = normalize_chromosome(parts[chromosome_idx]);
        let position: u64 = parts[position_idx]
            .parse()
            .with_context(|| format!("Invalid position: {}", parts[position_idx]))?;

        // Extract optional fields
        let rsid = column_mapping.get("rsid").and_then(|&idx| {
            if idx < parts.len() && parts[idx] != "." && !parts[idx].is_empty() {
                Some(parts[idx].to_string())
            } else {
                None
            }
        });

        let reference = column_mapping
            .get("reference")
            .and_then(|&idx| {
                if idx < parts.len() {
                    Some(parts[idx].to_string())
                } else {
                    None
                }
            })
            .unwrap_or_default();

        let alternate = column_mapping
            .get("alternate")
            .and_then(|&idx| {
                if idx < parts.len() {
                    Some(parts[idx].to_string())
                } else {
                    None
                }
            })
            .unwrap_or_default();

        let genotype = if let Some(&genotype_idx) = column_mapping.get("genotype") {
            if genotype_idx < parts.len() {
                parse_genotype(parts[genotype_idx], &reference, &[alternate.clone()])
            } else {
                Genotype::NoCall
            }
        } else {
            // Try to construct genotype from separate allele columns
            let allele1 = column_mapping.get("allele1").and_then(|&idx| {
                if idx < parts.len() {
                    Some(parts[idx])
                } else {
                    None
                }
            });
            let allele2 = column_mapping.get("allele2").and_then(|&idx| {
                if idx < parts.len() {
                    Some(parts[idx])
                } else {
                    None
                }
            });

            if let (Some(a1), Some(a2)) = (allele1, allele2) {
                parse_genotype(&format!("{}{}", a1, a2), &reference, &[alternate.clone()])
            } else {
                Genotype::NoCall
            }
        };

        let variant = Variant {
            chromosome,
            position,
            rsid,
            reference_allele: reference,
            alternate_alleles: if alternate.is_empty() {
                vec![]
            } else {
                vec![alternate]
            },
            genotype,
            quality: None,
            filter: None,
            info: Default::default(),
        };

        data.add_variant(variant);
        Ok(())
    }
}

impl super::GeneticDataParser for TsvParser {
    fn parse(&self, path: &Path) -> Result<ParsedGeneticData> {
        self.parse(path)
    }
}
