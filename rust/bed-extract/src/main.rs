use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::process;

// ── BED entry ──────────────────────────────────────────────────────────
struct BedEntry {
    start: usize,
    end: usize,
    strand: u8, // b'+' or b'-'
    period: usize,
}

// ── Bioinformatics helpers ─────────────────────────────────────────────

#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C',
        b'a' => b't', b't' => b'a', b'c' => b'g', b'g' => b'c',
        b'N' => b'N', b'n' => b'n',
        other => other,
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

fn gc_content(seq: &[u8]) -> f64 {
    if seq.is_empty() { return 0.0; }
    let gc = seq.iter().filter(|&&b| b == b'G' || b == b'C').count();
    gc as f64 / seq.len() as f64
}

fn shannon_entropy(seq: &[u8]) -> f64 {
    if seq.is_empty() { return 0.0; }
    let mut counts = [0u64; 4]; // A, C, G, T
    for &b in seq {
        match b {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            _ => {}
        }
    }
    let len = seq.len() as f64;
    let mut entropy = 0.0f64;
    for &c in &counts {
        if c > 0 {
            let p = c as f64 / len;
            entropy -= p * p.log2();
        }
    }
    entropy
}

fn calculate_pmatch(seq: &[u8], period: usize) -> u32 {
    if period == 0 || seq.len() < period { return 0; }
    let consensus = &seq[..period];
    let mut total_matches: usize = 0;
    let mut total_bases: usize = 0;

    for chunk in seq.chunks(period) {
        let cmp_len = chunk.len().min(period);
        for i in 0..cmp_len {
            if chunk[i] == consensus[i % period] {
                total_matches += 1;
            }
        }
        total_bases += cmp_len;
    }

    if total_bases == 0 { return 0; }
    ((100 * total_matches + total_bases / 2) / total_bases) as u32
}

// ── BED parsing ────────────────────────────────────────────────────────

fn load_bed(path: &str) -> io::Result<HashMap<String, Vec<BedEntry>>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut map: HashMap<String, Vec<BedEntry>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 { continue; }

        let chr = fields[0].split_whitespace().next().unwrap_or("").to_string();
        let start: usize = match fields[1].parse() { Ok(v) => v, Err(_) => continue };
        let end: usize = match fields[2].parse() { Ok(v) => v, Err(_) => continue };
        if start >= end { continue; }

        let strand = if fields.len() > 5 && fields[5] == "-" { b'-' } else { b'+' };
        let period: usize = if fields.len() > 3 { fields[3].parse().unwrap_or(1) } else { 1 };

        map.entry(chr).or_default().push(BedEntry { start, end, strand, period });
    }

    for entries in map.values_mut() {
        entries.sort_by_key(|e| e.start);
    }

    Ok(map)
}

// ── FASTA streaming + extraction ───────────────────────────────────────

fn extract_from_chr(
    name: &str,
    seq: &[u8],
    entries: &[BedEntry],
    sat_out: &mut BufWriter<File>,
    fasta_out: &mut Option<BufWriter<File>>,
    project: &str,
    trf_id: &mut u64,
) -> io::Result<usize> {
    let chr_len = seq.len();
    let mut count = 0usize;

    for entry in entries {
        if entry.end > chr_len { continue; }

        let raw = &seq[entry.start..entry.end];
        let extracted_seq: Vec<u8>;
        let upper: Vec<u8>;

        if entry.strand == b'-' {
            extracted_seq = reverse_complement(raw);
            upper = extracted_seq.iter().map(|&b| b.to_ascii_uppercase()).collect();
        } else {
            upper = raw.iter().map(|&b| b.to_ascii_uppercase()).collect();
        }

        let arr_len = upper.len();
        let period = entry.period;
        let n_copy = if period > 0 { format!("{:.1}", arr_len as f64 / period as f64) } else { "0".into() };
        let consensus_end = if period <= arr_len { period } else { arr_len };
        let consensus = &upper[..consensus_end];
        let arr_gc = gc_content(&upper);
        let con_gc = gc_content(consensus);
        let pmatch = calculate_pmatch(&upper, period);
        let pvar = 100u32.saturating_sub(pmatch);
        let entropy = shannon_entropy(&upper);

        *trf_id += 1;

        // SAT line: 18 fields
        write!(sat_out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t",
            project, trf_id, name, entry.start + 1, entry.end, period, n_copy, pmatch, pvar, entropy)?;
        sat_out.write_all(consensus)?;
        sat_out.write_all(b"\t")?;
        sat_out.write_all(&upper)?;
        write!(sat_out, "\t{:.2}\t{:.2}\t{}\t0\t\t\n", arr_gc, con_gc, arr_len)?;

        if let Some(ref mut fh) = fasta_out {
            write!(fh, ">{}_{}_{}_{}\n", name, entry.start, entry.end, period)?;
            fh.write_all(&upper)?;
            fh.write_all(b"\n")?;
        }

        count += 1;
    }
    Ok(count)
}

fn process_fasta<R: Read>(
    reader: R,
    bed: &HashMap<String, Vec<BedEntry>>,
    sat_out: &mut BufWriter<File>,
    fasta_out: &mut Option<BufWriter<File>>,
    project: &str,
) -> io::Result<usize> {
    let buf = BufReader::with_capacity(1 << 20, reader);
    let mut chr_name = String::new();
    let mut seq = Vec::<u8>::with_capacity(300_000_000); // pre-alloc ~300 MB
    let mut extracted = 0usize;
    let mut trf_id = 0u64;

    for line in buf.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !chr_name.is_empty() {
                if let Some(entries) = bed.get(&chr_name) {
                    let n = extract_from_chr(&chr_name, &seq, entries, sat_out, fasta_out, project, &mut trf_id)?;
                    extracted += n;
                }
            }
            chr_name = line[1..].split_whitespace().next().unwrap_or("").to_string();
            seq.clear();
            eprint!("\rProcessing {}...            ", chr_name);
        } else {
            seq.extend_from_slice(line.as_bytes());
        }
    }

    // Last chromosome
    if !chr_name.is_empty() {
        if let Some(entries) = bed.get(&chr_name) {
            let n = extract_from_chr(&chr_name, &seq, entries, sat_out, fasta_out, project, &mut trf_id)?;
            extracted += n;
        }
    }
    eprintln!();

    Ok(extracted)
}

// ── Main ───────────────────────────────────────────────────────────────

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 || args.len() > 6 {
        eprintln!("Usage: bed-extract <fasta[.gz]> <bed> <output.sat> [output.fasta] [project]");
        eprintln!("  Extracts sequences from FASTA at BED coordinates, outputs TRF/SAT format");
        process::exit(1);
    }

    let fasta_path = &args[1];
    let bed_path = &args[2];
    let sat_path = &args[3];
    let fasta_out_path = args.get(4).map(|s| s.as_str());
    let project = args.get(5).map(|s| s.as_str()).unwrap_or("FasTAN");

    // Load BED
    eprintln!("Loading BED file...");
    let bed = load_bed(bed_path).unwrap_or_else(|e| {
        eprintln!("Error reading BED: {}", e);
        process::exit(1);
    });
    let total_entries: usize = bed.values().map(|v| v.len()).sum();
    eprintln!("Loaded {} entries for {} chromosomes", total_entries, bed.len());

    // Open SAT output
    let sat_file = File::create(sat_path).unwrap_or_else(|e| {
        eprintln!("Error creating {}: {}", sat_path, e);
        process::exit(1);
    });
    let mut sat_out = BufWriter::with_capacity(1 << 20, sat_file);

    // SAT header
    writeln!(sat_out, "# bed-extract results in TRF format").unwrap();
    writeln!(sat_out, "# Source FASTA: {}", fasta_path).unwrap();
    writeln!(sat_out, "# Source BED: {}", bed_path).unwrap();
    writeln!(sat_out, "project\ttrf_id\ttrf_head\ttrf_l_ind\ttrf_r_ind\ttrf_period\ttrf_n_copy\ttrf_pmatch\ttrf_pvar\ttrf_entropy\ttrf_consensus\ttrf_array\ttrf_array_gc\ttrf_consensus_gc\ttrf_array_length\ttrf_joined\ttrf_family\ttrf_ref_annotation").unwrap();

    // Open optional FASTA output
    let mut fasta_out = fasta_out_path.map(|p| {
        let f = File::create(p).unwrap_or_else(|e| {
            eprintln!("Error creating {}: {}", p, e);
            process::exit(1);
        });
        BufWriter::with_capacity(1 << 20, f)
    });

    // Open FASTA input
    let fasta_file = File::open(fasta_path).unwrap_or_else(|e| {
        eprintln!("Error opening {}: {}", fasta_path, e);
        process::exit(1);
    });

    eprintln!("Processing FASTA...");
    let extracted = if fasta_path.ends_with(".gz") {
        let decoder = MultiGzDecoder::new(fasta_file);
        process_fasta(decoder, &bed, &mut sat_out, &mut fasta_out, project)
    } else {
        process_fasta(fasta_file, &bed, &mut sat_out, &mut fasta_out, project)
    }.unwrap_or_else(|e| {
        eprintln!("Error: {}", e);
        process::exit(1);
    });

    sat_out.flush().unwrap();
    if let Some(ref mut fh) = fasta_out {
        fh.flush().unwrap();
    }

    eprintln!("Extracted {} sequences", extracted);
}
