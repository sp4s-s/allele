use anyhow::{anyhow, Context, Result};
use std::path::Path;

use crate::parsers::{normalize_chromosome, open_file, parse_genotype, ParsedGeneticData};
use crate::types::*;

/// AncestryDNA parser for raw DNA data files
pub struct AncestryDNAParser;

impl AncestryDNAParser {
    pub fn new() -> Self {
        Self
    }

    pub fn parse(&self, path: &Path) -> Result<ParsedGeneticData> {
        let mut reader = open_file(path)?;
        let mut data = self.initialize_data(path)?;

        let mut line = String::new();
        while reader.read_line(&mut line)? > 0 {
            if line.starts_with('#') {
                // Skip comments
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
            FileFormat::AncestryDNA,
        ))
    }

    fn parse_data_line(&self, line: &str, data: &mut ParsedGeneticData) -> Result<()> {
        let line = line.trim();
        if line.is_empty() {
            return Ok(());
        }

        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 4 {
            return Err(anyhow!("Invalid AncestryDNA line: not enough columns"));
        }

        let rsid = parts[0].to_string();
        let chromosome = normalize_chromosome(parts[1]);
        let position: u64 = parts[2]
            .parse()
            .with_context(|| format!("Invalid position: {}", parts[2]))?;
        let allele1 = parts[3].chars().next().unwrap_or('N').to_string();
        let allele2 = parts[3].chars().nth(1).unwrap_or('N').to_string();
        let genotype_str = format!("{}{}", allele1, allele2);

        // For AncestryDNA files, we don't have explicit reference/alternate alleles
        let genotype = parse_genotype(&genotype_str, "", &[]);

        let variant = Variant {
            chromosome,
            position,
            rsid: Some(rsid),
            reference_allele: "".to_string(), // Not available in AncestryDNA format
            alternate_alleles: vec![],        // Not available in AncestryDNA format
            genotype,
            quality: None,
            filter: None,
            info: Default::default(),
        };

        data.add_variant(variant);
        Ok(())
    }
}

impl super::GeneticDataParser for AncestryDNAParser {
    fn parse(&self, path: &Path) -> Result<ParsedGeneticData> {
        self.parse(path)
    }
}
