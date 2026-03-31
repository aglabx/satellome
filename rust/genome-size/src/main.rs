use flate2::read::MultiGzDecoder;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;
use std::process;

struct SeqRecord {
    name: String,
    length: u64,
    offset: u64,      // byte offset of first base
    line_bases: u64,   // bases per line
    line_bytes: u64,   // bytes per line (bases + newline)
}

/// Parse FASTA, collect per-scaffold stats, write .fai and print summary.
/// For gzip input, offset/line_bases/line_bytes are approximate (not seekable).
fn process_fasta<R: Read>(reader: R, fai_path: &str) -> io::Result<()> {
    let buf = BufReader::with_capacity(1 << 20, reader); // 1 MB buffer
    let mut records: Vec<SeqRecord> = Vec::new();
    let mut current_name = String::new();
    let mut current_len: u64 = 0;
    let mut seq_offset: u64 = 0;
    let mut first_line_bases: u64 = 0;
    let mut first_line_bytes: u64 = 0;
    let mut byte_offset: u64 = 0;
    let mut seen_first_seq_line = false;

    for line_result in buf.lines() {
        let line = line_result?;
        let line_byte_len = line.len() as u64 + 1; // +1 for newline

        if line.starts_with('>') {
            if !current_name.is_empty() {
                records.push(SeqRecord {
                    name: current_name.clone(),
                    length: current_len,
                    offset: seq_offset,
                    line_bases: first_line_bases,
                    line_bytes: first_line_bytes,
                });
            }
            current_name = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            current_len = 0;
            seq_offset = byte_offset + line_byte_len;
            seen_first_seq_line = false;
        } else {
            current_len += line.len() as u64;
            if !seen_first_seq_line {
                first_line_bases = line.len() as u64;
                first_line_bytes = line_byte_len;
                seen_first_seq_line = true;
            }
        }
        byte_offset += line_byte_len;
    }

    if !current_name.is_empty() {
        records.push(SeqRecord {
            name: current_name,
            length: current_len,
            offset: seq_offset,
            line_bases: first_line_bases,
            line_bytes: first_line_bytes,
        });
    }

    // Write .fai file (for .gz, offsets refer to uncompressed stream)
    let mut fai = File::create(fai_path)?;
    for r in &records {
        writeln!(fai, "{}\t{}\t{}\t{}\t{}", r.name, r.length, r.offset, r.line_bases, r.line_bytes)?;
    }
    eprintln!("Wrote {}", fai_path);

    // Print per-scaffold sizes to stdout
    let mut total: u64 = 0;
    for r in &records {
        println!("{}\t{}", r.name, r.length);
        total += r.length;
    }
    println!("TOTAL\t{}\t{}", total, records.len());

    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: genome-size <file.(fa|fasta|fna)[.gz]>");
        eprintln!("Outputs: per-scaffold sizes (TSV) to stdout, creates .fai index");
        process::exit(1);
    }

    let path = &args[1];

    // Strip .gz to get the base FASTA path, then append .fai
    let base_path = if path.ends_with(".gz") {
        &path[..path.len() - 3]
    } else {
        path.as_str()
    };
    let fai_path = format!("{}.fai", base_path);
    if Path::new(&fai_path).exists() {
        eprintln!("{} already exists, reading from index", fai_path);
        // Parse existing .fai to print sizes
        let fai_file = File::open(&fai_path).unwrap_or_else(|e| {
            eprintln!("Error opening {}: {}", fai_path, e);
            process::exit(1);
        });
        let reader = BufReader::new(fai_file);
        let mut total: u64 = 0;
        let mut num_seq: u64 = 0;
        for line in reader.lines() {
            if let Ok(line) = line {
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() >= 2 {
                    let name = parts[0];
                    let len: u64 = parts[1].parse().unwrap_or(0);
                    println!("{}\t{}", name, len);
                    total += len;
                    num_seq += 1;
                }
            }
        }
        println!("TOTAL\t{}\t{}", total, num_seq);
        return;
    }

    let file = File::open(path).unwrap_or_else(|e| {
        eprintln!("Error opening {}: {}", path, e);
        process::exit(1);
    });

    let result = if path.ends_with(".gz") {
        let decoder = MultiGzDecoder::new(file);
        process_fasta(decoder, &fai_path)
    } else {
        process_fasta(file, &fai_path)
    };

    if let Err(e) = result {
        eprintln!("Error reading {}: {}", path, e);
        process::exit(1);
    }
}
