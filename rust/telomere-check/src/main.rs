use flate2::read::MultiGzDecoder;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write, BufWriter};
use std::process;

/// Default telomere motif for vertebrates: TTAGGG
/// At 5' (left) end we expect the reverse complement: CCCTAA repeats
/// At 3' (right) end we expect: TTAGGG repeats
/// This is because the telomere is double-stranded:
///   5'---TTAGGGTTAGGGTTAGGG-3'  (leading strand at right end)
///   3'---AATCCCAATCCCAATCCC-5'  (complement)
/// In the assembled sequence (always 5'->3'):
///   Left end:  CCCTAACCCTAACCCTAA...  (rev-comp of TTAGGG)
///   Right end: ...TTAGGGTTAGGGTTAGGG

const DEFAULT_MOTIF: &str = "TTAGGG";
const WINDOW_SIZE: usize = 10000; // bases to check at each end
const MIN_DENSITY: f64 = 0.15;   // minimum fraction of window with telomeric 4-mers to call PRESENT
const MIN_REPEATS: usize = 10;   // minimum total motif occurrences in window for PARTIAL

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C',
        b'a' => b't', b't' => b'a', b'c' => b'g', b'g' => b'c',
        other => other,
    }).collect()
}

/// Count maximum consecutive repeats of motif starting from the beginning of seq
fn count_forward_repeats(seq: &[u8], motif: &[u8]) -> usize {
    let mlen = motif.len();
    if mlen == 0 || seq.len() < mlen { return 0; }

    let mut count = 0;
    let mut pos = 0;

    // Skip leading N's
    while pos < seq.len() && (seq[pos] == b'N' || seq[pos] == b'n') {
        pos += 1;
    }

    // Allow 1 mismatch per repeat for degenerate telomeres
    while pos + mlen <= seq.len() {
        let chunk = &seq[pos..pos + mlen];
        if matches_with_mismatches(chunk, motif, 1) {
            count += 1;
            pos += mlen;
        } else {
            break;
        }
    }
    count
}

/// Count maximum consecutive repeats of motif ending at the end of seq
fn count_reverse_repeats(seq: &[u8], motif: &[u8]) -> usize {
    let mlen = motif.len();
    if mlen == 0 || seq.len() < mlen { return 0; }

    let mut count = 0;
    let mut end = seq.len();

    // Skip trailing N's
    while end > 0 && (seq[end - 1] == b'N' || seq[end - 1] == b'n') {
        end -= 1;
    }

    while end >= mlen {
        let chunk = &seq[end - mlen..end];
        if matches_with_mismatches(chunk, motif, 1) {
            count += 1;
            end -= mlen;
        } else {
            break;
        }
    }
    count
}

/// Check if two sequences match with at most max_mismatches differences
#[inline]
fn matches_with_mismatches(a: &[u8], b: &[u8], max_mismatches: usize) -> bool {
    if a.len() != b.len() { return false; }
    let mut mm = 0;
    for i in 0..a.len() {
        if a[i].to_ascii_uppercase() != b[i].to_ascii_uppercase() { mm += 1; }
        if mm > max_mismatches { return false; }
    }
    true
}

/// Count how many bases in window are covered by the telomere motif pattern.
/// Uses a sliding approach: for each position, check if the doubled motif
/// contains the local sequence as a substring (handles phase shifts and indels).
/// Returns (exact_motif_count, covered_bases).
fn telomere_coverage(seq: &[u8], motif: &[u8]) -> (usize, usize) {
    let mlen = motif.len();
    if mlen == 0 || seq.is_empty() { return (0, 0); }

    // Build doubled motif for phase-independent matching
    let mut doubled = Vec::with_capacity(mlen * 2);
    doubled.extend_from_slice(motif);
    doubled.extend_from_slice(motif);

    // Count exact motif hits
    let mut exact = 0;
    for i in 0..=seq.len().saturating_sub(mlen) {
        if matches_with_mismatches(&seq[i..i + mlen], motif, 1) {
            exact += 1;
        }
    }

    // Count bases covered: check 4-mers from the motif cycle
    // A base is "telomeric" if a 4-mer starting at that position
    // appears in the doubled motif
    let kmer_len = 4usize.min(mlen);
    let mut covered = 0;
    for i in 0..=seq.len().saturating_sub(kmer_len) {
        let kmer = &seq[i..i + kmer_len];
        // Check if this kmer appears anywhere in the doubled motif
        let mut found = false;
        for j in 0..=doubled.len().saturating_sub(kmer_len) {
            if kmer.eq_ignore_ascii_case(&doubled[j..j + kmer_len]) {
                found = true;
                break;
            }
        }
        if found {
            covered += 1;
        }
    }

    (exact, covered)
}

/// Find interstitial telomere sites (ITS) - telomere motifs away from chromosome ends
/// Returns list of (start, end, count) for regions with >= min_its consecutive motifs
fn find_interstitial_sites(seq: &[u8], motif: &[u8], revcomp: &[u8], min_its: usize) -> Vec<(usize, usize, usize)> {
    let mlen = motif.len();
    if mlen == 0 || seq.len() < mlen { return vec![]; }

    // Skip the terminal windows
    let skip = WINDOW_SIZE.min(seq.len() / 2);
    if seq.len() <= skip * 2 { return vec![]; }

    let interior = &seq[skip..seq.len() - skip];
    let mut sites = Vec::new();
    let mut i = 0;
    while i + mlen <= interior.len() {
        let chunk = &interior[i..i + mlen];
        if matches_with_mismatches(chunk, motif, 1) || matches_with_mismatches(chunk, revcomp, 1) {
            let start = i;
            let mut count = 0;
            while i + mlen <= interior.len() {
                let c = &interior[i..i + mlen];
                if matches_with_mismatches(c, motif, 1) || matches_with_mismatches(c, revcomp, 1) {
                    count += 1;
                    i += mlen;
                } else {
                    break;
                }
            }
            if count >= min_its {
                sites.push((start + skip, start + skip + count * mlen, count));
            }
        } else {
            i += 1;
        }
    }
    sites
}

struct TelomereResult {
    left_status: String,
    right_status: String,
    t2t: bool,
}

fn check_telomeres(seq: &[u8], motif: &[u8]) -> (usize, usize, f64, String, usize, usize, f64, String, bool) {
    let revcomp_motif = reverse_complement(motif);
    let window = WINDOW_SIZE.min(seq.len());

    let left_window = &seq[..window];
    let right_window = &seq[seq.len() - window..];

    // Left end: expect rev-comp of motif (CCCTAA for vertebrates)
    let left_consecutive = count_forward_repeats(left_window, &revcomp_motif);
    let (left_exact, left_covered) = telomere_coverage(left_window, &revcomp_motif);
    let left_density = left_covered as f64 / window as f64;

    // Right end: expect forward motif (TTAGGG for vertebrates)
    let right_consecutive = count_reverse_repeats(right_window, motif);
    let (right_exact, right_covered) = telomere_coverage(right_window, motif);
    let right_density = right_covered as f64 / window as f64;

    let left_status = if left_density >= MIN_DENSITY {
        "PRESENT".to_string()
    } else if left_exact >= MIN_REPEATS {
        "PARTIAL".to_string()
    } else {
        "ABSENT".to_string()
    };

    let right_status = if right_density >= MIN_DENSITY {
        "PRESENT".to_string()
    } else if right_exact >= MIN_REPEATS {
        "PARTIAL".to_string()
    } else {
        "ABSENT".to_string()
    };

    let t2t = left_status == "PRESENT" && right_status == "PRESENT";

    (left_consecutive, left_exact, left_density, left_status,
     right_consecutive, right_exact, right_density, right_status, t2t)
}

fn process_fasta<R: Read>(reader: R, motif: &[u8], out: &mut BufWriter<File>,
                          its_out: &mut BufWriter<File>) -> io::Result<Vec<TelomereResult>> {
    let buf = BufReader::with_capacity(1 << 20, reader);
    let revcomp = reverse_complement(motif);
    let mut chr_name = String::new();
    let mut seq = Vec::<u8>::with_capacity(300_000_000);
    let mut results = Vec::new();

    for line in buf.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !chr_name.is_empty() {
                let upper: Vec<u8> = seq.iter().map(|&b| b.to_ascii_uppercase()).collect();
                let (lc, lt, ld, ls, rc, rt, rd, rs, t2t) = check_telomeres(&upper, motif);

                writeln!(out, "{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{:.3}\t{}\t{}",
                    chr_name, upper.len(), lc, lt, ld, ls, rc, rt, rd, rs,
                    if t2t { "T2T" } else { "incomplete" })?;

                // Find interstitial telomere sites
                let its = find_interstitial_sites(&upper, motif, &revcomp, MIN_REPEATS as usize);
                for (start, end, count) in &its {
                    writeln!(its_out, "{}\t{}\t{}\t{}", chr_name, start, end, count)?;
                }

                results.push(TelomereResult { left_status: ls, right_status: rs, t2t });
            }
            chr_name = line[1..].split_whitespace().next().unwrap_or("").to_string();
            seq.clear();
            eprint!("\r{}            ", chr_name);
        } else {
            seq.extend_from_slice(line.as_bytes());
        }
    }

    // Last chromosome
    if !chr_name.is_empty() {
        let upper: Vec<u8> = seq.iter().map(|&b| b.to_ascii_uppercase()).collect();
        let (lc, lt, ld, ls, rc, rt, rd, rs, t2t) = check_telomeres(&upper, motif);
        writeln!(out, "{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{:.3}\t{}\t{}",
            chr_name, upper.len(), lc, lt, ld, ls, rc, rt, rd, rs,
            if t2t { "T2T" } else { "incomplete" })?;
        let its = find_interstitial_sites(&upper, motif, &revcomp, MIN_REPEATS as usize);
        for (start, end, count) in &its {
            writeln!(its_out, "{}\t{}\t{}\t{}", chr_name, start, end, count)?;
        }
        results.push(TelomereResult { left_status: ls, right_status: rs, t2t });
    }
    eprintln!();

    Ok(results)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 || args.len() > 4 {
        eprintln!("Usage: telomere-check <fasta[.gz]> <output.tsv> [motif]");
        eprintln!("");
        eprintln!("  motif: telomere repeat unit (default: TTAGGG for vertebrates)");
        eprintln!("         Common motifs:");
        eprintln!("           TTAGGG   - vertebrates, fungi, many metazoans");
        eprintln!("           TTTAGGG  - plants (Arabidopsis)");
        eprintln!("           TTAGG    - insects");
        eprintln!("           TTAGGC   - nematodes (C. elegans)");
        eprintln!("           TTGGGG   - Tetrahymena");
        eprintln!("");
        eprintln!("  Checks both ends of each chromosome:");
        eprintln!("    Left (5'):  expects reverse complement of motif");
        eprintln!("    Right (3'): expects forward motif");
        eprintln!("  PRESENT: >{:.0}% of {}bp window covered by motif", MIN_DENSITY * 100.0, WINDOW_SIZE);
        eprintln!("  Also outputs interstitial telomere sites (ITS) to <output>.its.bed");
        process::exit(1);
    }

    let fasta_path = &args[1];
    let output_path = &args[2];
    let motif_str = args.get(3).map(|s| s.as_str()).unwrap_or(DEFAULT_MOTIF);
    let motif = motif_str.as_bytes();
    let revcomp = reverse_complement(motif);
    let revcomp_str = String::from_utf8_lossy(&revcomp);

    eprintln!("Motif: {} (left end expects: {})", motif_str, revcomp_str);

    let fasta_file = File::open(fasta_path).unwrap_or_else(|e| {
        eprintln!("Error opening {}: {}", fasta_path, e);
        process::exit(1);
    });

    let out_file = File::create(output_path).unwrap_or_else(|e| {
        eprintln!("Error creating {}: {}", output_path, e);
        process::exit(1);
    });
    let mut out = BufWriter::new(out_file);

    // ITS output file
    let its_path = output_path.replace(".tsv", ".its.bed");
    let its_file = File::create(&its_path).unwrap_or_else(|e| {
        eprintln!("Error creating {}: {}", its_path, e);
        process::exit(1);
    });
    let mut its_out = BufWriter::new(its_file);
    writeln!(its_out, "# Interstitial telomere sites (ITS)").unwrap();
    writeln!(its_out, "# chr\tstart\tend\trepeat_count").unwrap();

    // Header
    writeln!(out, "# telomere-check v0.1.0").unwrap();
    writeln!(out, "# motif: {} | rev-comp: {}", motif_str, revcomp_str).unwrap();
    writeln!(out, "# window: {} bp | min_density: {}", WINDOW_SIZE, MIN_DENSITY).unwrap();
    writeln!(out, "chr\tlength\tleft_consec\tleft_total\tleft_density\tleft_status\tright_consec\tright_total\tright_density\tright_status\tt2t_status").unwrap();

    let results = if fasta_path.ends_with(".gz") {
        let decoder = MultiGzDecoder::new(fasta_file);
        process_fasta(decoder, motif, &mut out, &mut its_out)
    } else {
        process_fasta(fasta_file, motif, &mut out, &mut its_out)
    }.unwrap_or_else(|e| {
        eprintln!("Error: {}", e);
        process::exit(1);
    });

    out.flush().unwrap();
    its_out.flush().unwrap();

    // Summary
    let total = results.len();
    let t2t_count = results.iter().filter(|r| r.t2t).count();
    let left_ok = results.iter().filter(|r| r.left_status == "PRESENT").count();
    let right_ok = results.iter().filter(|r| r.right_status == "PRESENT").count();

    eprintln!("Summary: {}/{} chromosomes T2T", t2t_count, total);
    eprintln!("  Left telomere present:  {}/{}", left_ok, total);
    eprintln!("  Right telomere present: {}/{}", right_ok, total);
    if t2t_count == total {
        eprintln!("  Assembly: FULLY T2T");
    } else {
        eprintln!("  Assembly: INCOMPLETE ({} chromosomes missing telomeres)", total - t2t_count);
    }
}
