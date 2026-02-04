use anyhow::{Context, Result};
use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

use crate::parsers::FileParser;
use crate::types::FileFormat;

/// File discovery system for finding genetic data files
pub struct FileDiscovery {
    recursive: bool,
}

impl FileDiscovery {
    pub fn new(recursive: bool) -> Self {
        Self { recursive }
    }

    /// Discover genetic data files based on patient file and comparison paths
    pub fn discover(&self, patient_path: &Path, compare_paths: &[PathBuf]) -> Result<Vec<PathBuf>> {
        let mut files = Vec::new();

        // Add patient file if it exists
        if patient_path.exists() {
            files.push(patient_path.to_path_buf());
        }

        // Process comparison paths
        for path in compare_paths {
            if path.is_file() {
                files.push(path.clone());
            } else if path.is_dir() {
                let dir_files = self.discover_in_directory(path)?;
                files.extend(dir_files);
            }
        }

        // If no files found and recursive mode, search current directory
        if files.is_empty() && self.recursive {
            let current_dir = std::env::current_dir().context("Failed to get current directory")?;
            let dir_files = self.discover_in_directory(&current_dir)?;
            files.extend(dir_files);
        }

        // Remove duplicates while preserving order
        let mut seen = HashSet::new();
        files.retain(|path| seen.insert(path.clone()));

        // Filter to only include actual genetic data files
        self.filter_genetic_files(files)
    }

    fn discover_in_directory(&self, dir: &Path) -> Result<Vec<PathBuf>> {
        let mut files = Vec::new();

        if self.recursive {
            // Recursive search
            for entry in WalkDir::new(dir)
                .follow_links(true)
                .into_iter()
                .filter_map(|e| e.ok())
            {
                let path = entry.path();
                if path.is_file() {
                    // Check if file might be a genetic data file
                    if self.is_potential_genetic_file(path) {
                        files.push(path.to_path_buf());
                    }
                }
            }
        } else {
            // Non-recursive search
            let entries = fs::read_dir(dir)
                .with_context(|| format!("Failed to read directory: {}", dir.display()))?;

            for entry in entries {
                let entry = entry.with_context(|| {
                    format!("Failed to read directory entry in: {}", dir.display())
                })?;
                let path = entry.path();

                if path.is_file() && self.is_potential_genetic_file(&path) {
                    files.push(path);
                }
            }
        }

        Ok(files)
    }

    fn is_potential_genetic_file(&self, path: &Path) -> bool {
        if let Some(ext) = path.extension() {
            let ext_str = ext.to_string_lossy().to_lowercase();

            // Check for common genetic data file extensions
            match ext_str.as_str() {
                "vcf" | "bcf" | "bed" | "fa" | "fasta" | "fna" | "fq" | "fastq" | "txt" | "csv"
                | "tsv" | "ped" | "map" | "bim" | "fam" | "gvf" | "gff" | "gff3" | "sam"
                | "bam" | "cram" | "gb" | "gbk" | "gz" | "bz2" | "xz" => true,
                _ => {
                    // For unknown extensions, check if it's a text file that might contain genetic data
                    self.is_text_file_with_genetic_patterns(path)
                }
            }
        } else {
            // No extension - check content
            self.is_text_file_with_genetic_patterns(path)
        }
    }

    fn is_text_file_with_genetic_patterns(&self, path: &Path) -> bool {
        // Try to read first few lines to check for genetic data patterns
        if let Ok(file) = std::fs::File::open(path) {
            use std::io::{BufRead, BufReader};
            let reader = BufReader::new(file);
            let mut lines = reader.lines();

            // Check first 5 lines for genetic data patterns
            for _ in 0..5 {
                if let Some(Ok(line)) = lines.next() {
                    // Check for common genetic data headers/patterns
                    let line_lower = line.to_lowercase();

                    if line_lower.contains("rsid")
                        && (line_lower.contains("chromosome") || line_lower.contains("chr"))
                        && line_lower.contains("position")
                    {
                        return true;
                    }

                    if line.starts_with("##fileformat=VCF") {
                        return true;
                    }

                    if line.starts_with("#") && line.contains("AncestryDNA") {
                        return true;
                    }

                    if line.contains("# This data file generated by 23andMe") {
                        return true;
                    }
                }
            }
        }

        false
    }

    fn filter_genetic_files(&self, files: Vec<PathBuf>) -> Result<Vec<PathBuf>> {
        let parser = FileParser::new();
        let mut valid_files = Vec::new();

        for file in files {
            // Quick check based on extension
            if let Some(ext) = file.extension() {
                let format = FileFormat::from_extension(&ext.to_string_lossy());
                if format.is_genetic_format() {
                    valid_files.push(file);
                    continue;
                }
            }

            // For unknown formats, try to parse to confirm it's genetic data
            match parser.parse(&file) {
                Ok(_) => valid_files.push(file),
                Err(_) => {
                    // Not a valid genetic data file, skip it
                    continue;
                }
            }
        }

        Ok(valid_files)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::TempDir;

    #[test]
    fn test_file_discovery() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let dir_path = temp_dir.path();

        // Create test files
        let vcf_path = dir_path.join("test.vcf");
        let mut vcf_file = File::create(&vcf_path)?;
        writeln!(vcf_file, "##fileformat=VCFv4.2")?;
        writeln!(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
        writeln!(vcf_file, "1\t100\t.\tA\tT\t30\tPASS\t.")?;

        let txt_path = dir_path.join("data.txt");
        let mut txt_file = File::create(&txt_path)?;
        writeln!(txt_file, "rsid\tchromosome\tposition\tgenotype")?;
        writeln!(txt_file, "rs123\t1\t1000\tAA")?;

        let invalid_path = dir_path.join("invalid.txt");
        let mut invalid_file = File::create(&invalid_path)?;
        writeln!(invalid_file, "This is not genetic data")?;

        // Test discovery
        let discovery = FileDiscovery::new(false);
        let files = discovery.discover_in_directory(dir_path)?;

        // Should find 2 valid genetic files
        assert_eq!(files.len(), 2);
        assert!(files.contains(&vcf_path));
        assert!(files.contains(&txt_path));
        assert!(!files.contains(&invalid_path));

        Ok(())
    }
}
