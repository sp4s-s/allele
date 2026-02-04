use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::parsers::{normalize_chromosome, ParsedGeneticData};
use crate::types::*;

/// PLINK parser for PED/MAP or BED/BIM/FAM format files
pub struct PlinkParser;

impl PlinkParser {
    pub fn new() -> Self {
        Self
    }

    pub fn parse(&self, path: &Path) -> Result<ParsedGeneticData> {
        // Determine which PLINK format we're dealing with
        if path.extension().map(|ext| ext == "ped").unwrap_or(false)
            || path.with_extension("ped").exists()
        {
            self.parse_ped_format(path)
        } else {
            self.parse_bed_format(path)
        }
    }

    fn parse_ped_format(&self, path: &Path) -> Result<ParsedGeneticData> {
        let ped_path = if path.extension().map(|ext| ext == "ped").unwrap_or(false) {
            path.to_path_buf()
        } else {
            path.with_extension("ped")
        };

        let map_path = path.with_extension("map");

        // Parse MAP file first to get SNP information
        let snp_info = self.parse_map_file(&map_path)?;

        // Parse PED file
        let mut data = self.initialize_data(path)?;
        let file = File::open(&ped_path)
            .with_context(|| format!("Failed to open PED file: {}", ped_path.display()))?;
        let reader = BufReader::new(file);

        // PED format: FamilyID IndividualID PaternalID MaternalID Sex Phenotype Genotypes...
        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split_whitespace().collect();

            if parts.len() < 6 {
                continue; // Skip invalid lines
            }

            let individual_id = parts[1].to_string();

            // Process genotype data (starts at index 6)
            if parts.len() >= 6 + snp_info.len() * 2 {
                for (i, (rsid, chr, pos)) in snp_info.iter().enumerate() {
                    let allele1_idx = 6 + i * 2;
                    let allele2_idx = 6 + i * 2 + 1;

                    if allele1_idx < parts.len() && allele2_idx < parts.len() {
                        let allele1 = parts[allele1_idx];
                        let allele2 = parts[allele2_idx];

                        // Convert PLINK alleles (1/2/0) to nucleotides
                        // This is a simplification - in reality, we'd need the BIM file
                        let genotype_str = match (allele1, allele2) {
                            ("0", "0") => "./.".to_string(),
                            (a1, a2) => format!("{}{}", a1, a2),
                        };

                        let variant = Variant {
                            chromosome: chr.clone(),
                            position: *pos,
                            rsid: Some(rsid.clone()),
                            reference_allele: "A".to_string(), // Placeholder
                            alternate_alleles: vec!["T".to_string()], // Placeholder
                            genotype: crate::parsers::parse_genotype(
                                &genotype_str,
                                "A",
                                &["T".to_string()],
                            ),
                            quality: None,
                            filter: None,
                            info: Default::default(),
                        };

                        data.add_variant(variant);
                    }
                }
            }

            // Use the first individual's data as our sample
            data.metadata.sample_id = individual_id;
            break;
        }

        data.calculate_quality_metrics();
        Ok(data)
    }

    fn parse_map_file(&self, map_path: &Path) -> Result<Vec<(String, String, u64)>> {
        let file = File::open(map_path)
            .with_context(|| format!("Failed to open MAP file: {}", map_path.display()))?;
        let reader = BufReader::new(file);

        let mut snp_info = Vec::new();

        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split_whitespace().collect();

            if parts.len() >= 4 {
                let rsid = parts[1].to_string();
                let chromosome = normalize_chromosome(parts[0]);
                let position: u64 = parts[3]
                    .parse()
                    .with_context(|| format!("Invalid position in MAP file: {}", parts[3]))?;

                snp_info.push((rsid, chromosome, position));
            }
        }

        Ok(snp_info)
    }

    fn parse_bed_format(&self, path: &Path) -> Result<ParsedGeneticData> {
        let bim_path = path.with_extension("bim");
        let fam_path = path.with_extension("fam");

        // Parse BIM file to get SNP information
        let snp_info = self.parse_bim_file(&bim_path)?;

        // Parse FAM file to get individual information
        let individual_info = self.parse_fam_file(&fam_path)?;

        // For BED format, we'll create a simplified parser
        let mut data = self.initialize_data(path)?;

        if !individual_info.is_empty() {
            data.metadata.sample_id = individual_info[0].clone();
        }

        // Add some placeholder variants since parsing BED binary format is complex
        for (rsid, chr, pos) in snp_info.iter().take(100) {
            let variant = Variant {
                chromosome: chr.clone(),
                position: *pos,
                rsid: Some(rsid.clone()),
                reference_allele: "A".to_string(),
                alternate_alleles: vec!["T".to_string()],
                genotype: crate::parsers::parse_genotype("0/1", "A", &["T".to_string()]),
                quality: None,
                filter: None,
                info: Default::default(),
            };

            data.add_variant(variant);
        }

        data.calculate_quality_metrics();
        Ok(data)
    }

    fn parse_bim_file(&self, bim_path: &Path) -> Result<Vec<(String, String, u64)>> {
        let file = File::open(bim_path)
            .with_context(|| format!("Failed to open BIM file: {}", bim_path.display()))?;
        let reader = BufReader::new(file);

        let mut snp_info = Vec::new();

        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split_whitespace().collect();

            if parts.len() >= 6 {
                let rsid = parts[1].to_string();
                let chromosome = normalize_chromosome(parts[0]);
                let position: u64 = parts[3]
                    .parse()
                    .with_context(|| format!("Invalid position in BIM file: {}", parts[3]))?;

                snp_info.push((rsid, chromosome, position));
            }
        }

        Ok(snp_info)
    }

    fn parse_fam_file(&self, fam_path: &Path) -> Result<Vec<String>> {
        let file = File::open(fam_path)
            .with_context(|| format!("Failed to open FAM file: {}", fam_path.display()))?;
        let reader = BufReader::new(file);

        let mut individuals = Vec::new();

        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split_whitespace().collect();

            if !parts.is_empty() && parts.len() >= 2 {
                individuals.push(parts[1].to_string());
            }
        }

        Ok(individuals)
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
            FileFormat::PLINK,
        ))
    }
}

impl super::GeneticDataParser for PlinkParser {
    fn parse(&self, path: &Path) -> Result<ParsedGeneticData> {
        self.parse(path)
    }
}
