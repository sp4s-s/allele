use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::parsers::{normalize_chromosome, open_file, parse_genotype, ParsedGeneticData};
use crate::types::*;

/// VCF parser for Variant Call Format files
pub struct VcfParser;

impl VcfParser {
    pub fn new() -> Self {
        Self
    }

    pub fn parse(&self, path: &Path) -> Result<ParsedGeneticData> {
        let mut reader = open_file(path)?;
        let mut data = self.initialize_data(path)?;

        let mut line = String::new();
        while reader.read_line(&mut line)? > 0 {
            if line.starts_with('#') {
                self.parse_header_line(&line, &mut data)?;
                line.clear();
                continue;
            }

            self.parse_variant_line(&line, &mut data)?;
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
            FileFormat::VCF,
        ))
    }

    fn parse_header_line(&self, line: &str, data: &mut ParsedGeneticData) -> Result<()> {
        if line.starts_with("##fileformat=") {
            // Already handled in format detection
        } else if line.starts_with("##reference=") {
            let reference = line.strip_prefix("##reference=").unwrap_or("").trim();
            data.metadata.genome_build = reference.to_string();
        } else if line.starts_with("##FORMAT=<ID=GT,") {
            // GT format field definition
        } else if line.starts_with("#CHROM") {
            // Column header line - extract sample IDs if present
            let parts: Vec<&str> = line.trim_start_matches('#').split('\t').collect();
            if parts.len() >= 10 {
                data.metadata.sample_id = parts[9].to_string();
            }
        }
        Ok(())
    }

    fn parse_variant_line(&self, line: &str, data: &mut ParsedGeneticData) -> Result<()> {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 8 {
            return Err(anyhow!("Invalid VCF line: not enough columns"));
        }

        let chromosome = normalize_chromosome(parts[0]);
        let position: u64 = parts[1]
            .parse()
            .with_context(|| format!("Invalid position: {}", parts[1]))?;
        let rsid = if parts[2] == "." {
            None
        } else {
            Some(parts[2].to_string())
        };
        let reference = parts[3].to_string();
        let alternate = parts[4].to_string();
        let quality: Option<f32> = if parts[5] == "." {
            None
        } else {
            Some(parts[5].parse().unwrap_or(0.0))
        };
        let filter = if parts[6] == "." {
            None
        } else {
            Some(parts[6].to_string())
        };

        // Parse INFO field
        let info_map = self.parse_info_field(parts[7]);

        // Parse genotype if FORMAT and sample columns exist
        let mut genotype = Genotype::NoCall;
        if parts.len() >= 10 {
            genotype = self.parse_genotype_field(parts[8], parts[9], &reference, &alternate)?;
        }

        let variant = Variant {
            chromosome,
            position,
            rsid,
            reference_allele: reference,
            alternate_alleles: if alternate == "." {
                vec![]
            } else {
                alternate.split(',').map(|s| s.to_string()).collect()
            },
            genotype,
            quality,
            filter,
            info: info_map,
        };

        data.add_variant(variant);
        Ok(())
    }

    fn parse_info_field(&self, info_str: &str) -> HashMap<String, String> {
        let mut info_map = HashMap::new();
        if info_str == "." {
            return info_map;
        }

        for entry in info_str.split(';') {
            if let Some(eq_pos) = entry.find('=') {
                let key = entry[..eq_pos].to_string();
                let value = entry[eq_pos + 1..].to_string();
                info_map.insert(key, value);
            } else {
                info_map.insert(entry.to_string(), "true".to_string());
            }
        }

        info_map
    }

    fn parse_genotype_field(
        &self,
        format_str: &str,
        sample_str: &str,
        reference: &str,
        alternate: &str,
    ) -> Result<Genotype> {
        let format_fields: Vec<&str> = format_str.split(':').collect();
        let sample_fields: Vec<&str> = sample_str.split(':').collect();

        // Find GT field index
        let gt_index = format_fields.iter().position(|&f| f == "GT");
        if gt_index.is_none() {
            return Ok(Genotype::NoCall);
        }

        let gt_index = gt_index.unwrap();
        if gt_index >= sample_fields.len() {
            return Ok(Genotype::NoCall);
        }

        let gt_value = sample_fields[gt_index];
        Ok(parse_genotype(
            gt_value,
            reference,
            &[alternate.to_string()],
        ))
    }
}

impl super::GeneticDataParser for VcfParser {
    fn parse(&self, path: &Path) -> Result<ParsedGeneticData> {
        self.parse(path)
    }
}
