use anyhow::{anyhow, Context, Result};
use std::path::Path;

use crate::parsers::{normalize_chromosome, open_file, parse_genotype, ParsedGeneticData};
use crate::types::*;

/// MyHeritage parser for MyHeritage DNA raw data files
pub struct MyHeritageParser;

impl MyHeritageParser {
    pub fn new() -> Self {
        Self
    }

    pub fn parse(&self, path: &Path) -> Result<ParsedGeneticData> {
        let mut reader = open_file(path)?;
        let mut data = self.initialize_data(path)?;

        let mut line = String::new();
        while reader.read_line(&mut line)? > 0 {
            if line.starts_with('#') || line.starts_with("RSID") {
                // Skip comments and header
                line.clear();
                continue;
            }

            self.parse_data_line(&line, &mut data)?;
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

        Ok(ParsedGeneticData::new(
            filename.clone(),
            path.to_string_lossy().to_string(),
            FileFormat::MyHeritage,
        ))
    }

    fn parse_data_line(&self, line: &str, data: &mut ParsedGeneticData) -> Result<()> {
        let line = line.trim();
        if line.is_empty() {
            return Ok(());
        }

        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 4 {
            return Err(anyhow!("Invalid MyHeritage line: not enough columns"));
        }

        let rsid = parts[0].trim_matches('"').to_string();
        let chromosome = normalize_chromosome(parts[1].trim_matches('"'));
        let position: u64 = parts[2]
            .trim_matches('"')
            .parse()
            .with_context(|| format!("Invalid position: {}", parts[2]))?;
        let result = parts[3].trim_matches('"');

        // For MyHeritage files, we don't have explicit reference/alternate alleles
        let genotype = parse_genotype(result, "", &[]);

        let variant = Variant {
            chromosome,
            position,
            rsid: Some(rsid),
            reference_allele: "".to_string(), // Not available in MyHeritage format
            alternate_alleles: vec![],        // Not available in MyHeritage format
            genotype,
            quality: None,
            filter: None,
            info: Default::default(),
        };

        data.add_variant(variant);
        Ok(())
    }
}

impl super::GeneticDataParser for MyHeritageParser {
    fn parse(&self, path: &Path) -> Result<ParsedGeneticData> {
        self.parse(path)
    }
}
