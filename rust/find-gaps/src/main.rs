use flate2::read::MultiGzDecoder;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write, BufWriter};
use std::process;

/// Stream FASTA, find runs of 'N'/'n', output BED-like: chr\tstart\tend\tlength
/// Skips scaffolds shorter than min_length.
fn find_gaps<R: Read>(reader: R, min_length: usize, out: &mut BufWriter<File>) -> io::Result<usize> {
    let buf = BufReader::with_capacity(1 << 20, reader);
    let mut chr_name = String::new();
    let mut pos: usize = 0;
    let mut gap_start: Option<usize> = None;
    let mut total_gaps = 0usize;
    let mut chr_len: usize = 0;
    let mut pending_gaps: Vec<(usize, usize)> = Vec::new();

    for line in buf.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !chr_name.is_empty() {
                if let Some(gs) = gap_start {
                    pending_gaps.push((gs, pos));
                    gap_start = None;
                }
                if chr_len >= min_length {
                    for &(start, end) in &pending_gaps {
                        writeln!(out, "{}\t{}\t{}\t{}", chr_name, start, end, end - start)?;
                        total_gaps += 1;
                    }
                }
                pending_gaps.clear();
            }
            chr_name = line[1..].split_whitespace().next().unwrap_or("").to_string();
            pos = 0;
            chr_len = 0;
            eprint!("\r{}            ", chr_name);
        } else {
            let bytes = line.as_bytes();
            chr_len += bytes.len();

            for &b in bytes {
                if b == b'N' || b == b'n' {
                    if gap_start.is_none() {
                        gap_start = Some(pos);
                    }
                } else if let Some(gs) = gap_start {
                    pending_gaps.push((gs, pos));
                    gap_start = None;
                }
                pos += 1;
            }
        }
    }

    // Last chromosome
    if !chr_name.is_empty() {
        if let Some(gs) = gap_start {
            pending_gaps.push((gs, pos));
        }
        if chr_len >= min_length {
            for &(start, end) in &pending_gaps {
                writeln!(out, "{}\t{}\t{}\t{}", chr_name, start, end, end - start)?;
                total_gaps += 1;
            }
        }
    }
    eprintln!();

    Ok(total_gaps)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 || args.len() > 4 {
        eprintln!("Usage: find-gaps <fasta[.gz]> <output.bed> [min_scaffold_length]");
        eprintln!("  Finds all N-gaps in FASTA, outputs BED: chr\\tstart\\tend\\tlength");
        process::exit(1);
    }

    let fasta_path = &args[1];
    let output_path = &args[2];
    let min_length: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(100000);

    let fasta_file = File::open(fasta_path).unwrap_or_else(|e| {
        eprintln!("Error opening {}: {}", fasta_path, e);
        process::exit(1);
    });

    let out_file = File::create(output_path).unwrap_or_else(|e| {
        eprintln!("Error creating {}: {}", output_path, e);
        process::exit(1);
    });
    let mut out = BufWriter::new(out_file);

    eprintln!("Finding gaps (min scaffold: {} bp)...", min_length);
    let total = if fasta_path.ends_with(".gz") {
        let decoder = MultiGzDecoder::new(fasta_file);
        find_gaps(decoder, min_length, &mut out)
    } else {
        find_gaps(fasta_file, min_length, &mut out)
    }.unwrap_or_else(|e| {
        eprintln!("Error: {}", e);
        process::exit(1);
    });

    out.flush().unwrap();
    eprintln!("Found {} gaps", total);
}
