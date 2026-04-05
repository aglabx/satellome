use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write, BufWriter};
use std::process;

// ── DNA helpers ────────────────────────────────────────────────────────

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C',
        b'a' => b't', b't' => b'a', b'c' => b'g', b'g' => b'c',
        other => other,
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Generate all cyclic rotations of a sequence
fn rotations(seq: &[u8]) -> Vec<Vec<u8>> {
    let n = seq.len();
    if n == 0 { return vec![]; }
    let doubled: Vec<u8> = seq.iter().chain(seq.iter()).copied().collect();
    (0..n).map(|i| doubled[i..i+n].to_vec()).collect()
}

/// Canonical form: lexicographically smallest among all rotations of
/// the sequence and its reverse complement. This normalizes for:
/// - Cyclic rotations (ATCG == TCGA == CGAT == GATC)
/// - Strand (ATCG == CGAT on opposite strand)
fn canonical_form(seq: &[u8]) -> Vec<u8> {
    let upper: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
    let rc = reverse_complement(&upper);

    let mut best = upper.clone();
    for rot in rotations(&upper) {
        if rot < best { best = rot; }
    }
    for rot in rotations(&rc) {
        if rot < best { best = rot; }
    }
    best
}

/// Check if seq1 is a multiple of seq2 (or vice versa)
/// e.g., ATCGATCG (period 8) is 2x of ATCG (period 4)
/// Returns the shorter canonical period if one divides the other
fn find_base_period(seq: &[u8]) -> Vec<u8> {
    let n = seq.len();
    for divisor in 1..=n/2 {
        if n % divisor != 0 { continue; }
        let unit = &seq[..divisor];
        let repeats = n / divisor;
        let mut is_multiple = true;
        for r in 1..repeats {
            let chunk = &seq[r*divisor..(r+1)*divisor];
            if chunk != unit {
                is_multiple = false;
                break;
            }
        }
        if is_multiple {
            return unit.to_vec();
        }
    }
    seq.to_vec()
}

/// Edit distance between two sequences (Levenshtein)
fn edit_distance(a: &[u8], b: &[u8]) -> usize {
    let m = a.len();
    let n = b.len();
    let mut dp = vec![vec![0usize; n + 1]; m + 1];
    for i in 0..=m { dp[i][0] = i; }
    for j in 0..=n { dp[0][j] = j; }
    for i in 1..=m {
        for j in 1..=n {
            let cost = if a[i-1] == b[j-1] { 0 } else { 1 };
            dp[i][j] = (dp[i-1][j] + 1)
                .min(dp[i][j-1] + 1)
                .min(dp[i-1][j-1] + cost);
        }
    }
    dp[m][n]
}

/// Minimum edit distance considering all rotations and rev-comp
fn min_rotation_distance(a: &[u8], b: &[u8]) -> usize {
    if a.len() != b.len() {
        // Different lengths — try if one is a multiple of the other
        return edit_distance(a, b);
    }
    let rc_b = reverse_complement(b);
    let mut best = edit_distance(a, b);
    for rot in rotations(b) {
        let d = edit_distance(a, &rot);
        if d < best { best = d; }
        if best == 0 { return 0; }
    }
    for rot in rotations(&rc_b) {
        let d = edit_distance(a, &rot);
        if d < best { best = d; }
        if best == 0 { return 0; }
    }
    best
}

// ── Array record ───────────────────────────────────────────────────────

struct ArrayRecord {
    array_id: String,
    chr: String,
    start: usize,
    end: usize,
    period: usize,
    consensus: Vec<u8>,      // cut_sequence from first monomer
    canonical: Vec<u8>,       // canonical form of base period
    base_period: Vec<u8>,     // shortest non-redundant unit
    family_id: Option<usize>,
}

// ── Main logic ─────────────────────────────────────────────────────────

fn load_arrays(monomers_path: &str, min_length: usize) -> io::Result<Vec<ArrayRecord>> {
    let file = File::open(monomers_path)?;
    let reader = BufReader::new(file);
    let mut arrays: HashMap<String, ArrayRecord> = HashMap::new();
    let mut header_map: HashMap<String, usize> = HashMap::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        if line_num == 0 {
            // Parse header
            for (i, field) in line.split('\t').enumerate() {
                header_map.insert(field.to_string(), i);
            }
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        let get = |name: &str| -> &str {
            header_map.get(name).and_then(|&i| fields.get(i)).copied().unwrap_or("")
        };

        let aid = get("array_id").to_string();
        let rtype = get("type");

        // Parse array_id: chr_start_end_length_period
        let parts: Vec<&str> = aid.split('_').collect();
        if parts.len() < 5 { continue; }
        let arr_len: usize = parts[3].parse().unwrap_or(0);
        if arr_len < min_length { continue; }

        if rtype == "monomer" && !arrays.contains_key(&aid) {
            let cut_seq = get("cut_sequence").as_bytes().to_vec();
            let upper: Vec<u8> = cut_seq.iter().map(|b| b.to_ascii_uppercase()).collect();
            if upper.is_empty() { continue; }
            let base = find_base_period(&upper);
            let canon = canonical_form(&base);

            let chr = parts[0].to_string();
            let start = parts[1].parse().unwrap_or(0);
            let end = parts[2].parse().unwrap_or(0);
            let period = parts[4].parse().unwrap_or(0);
            let aid_clone = aid.clone();

            arrays.insert(aid, ArrayRecord {
                array_id: aid_clone,
                chr, start, end, period,
                consensus: upper,
                base_period: base,
                canonical: canon,
                family_id: None,
            });
        }
    }

    let mut result: Vec<ArrayRecord> = arrays.into_values().collect();
    result.sort_by(|a, b| a.chr.cmp(&b.chr).then(a.start.cmp(&b.start)));
    Ok(result)
}

/// Cluster arrays into families by canonical monomer similarity
fn cluster_families(arrays: &mut [ArrayRecord], max_edit_dist: usize) {
    // Group by canonical form first (exact matches)
    let mut canon_groups: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();
    for (i, arr) in arrays.iter().enumerate() {
        canon_groups.entry(arr.canonical.clone()).or_default().push(i);
    }

    // Assign family IDs to exact canonical matches
    let mut family_id = 0usize;
    let mut canon_to_family: HashMap<Vec<u8>, usize> = HashMap::new();
    let canon_keys: Vec<Vec<u8>> = {
        let mut keys: Vec<_> = canon_groups.keys().cloned().collect();
        // Sort by group size descending (largest families first)
        keys.sort_by(|a, b| canon_groups[b].len().cmp(&canon_groups[a].len()));
        keys
    };

    for canon in &canon_keys {
        if canon_to_family.contains_key(canon) { continue; }

        canon_to_family.insert(canon.clone(), family_id);

        // Try to merge similar canonical forms
        for other_canon in &canon_keys {
            if canon_to_family.contains_key(other_canon) { continue; }
            if other_canon == canon { continue; }

            // Check if base periods are compatible (same length or multiple)
            let a_bp = &arrays[canon_groups[canon][0]].base_period;
            let b_bp = &arrays[canon_groups[other_canon][0]].base_period;

            // Only compare if similar period length (within 2x)
            let len_ratio = a_bp.len().max(b_bp.len()) as f64 / a_bp.len().min(b_bp.len()).max(1) as f64;
            if len_ratio > 2.0 { continue; }

            // Compare using rotation-aware edit distance
            let dist = min_rotation_distance(a_bp, b_bp);
            let max_len = a_bp.len().max(b_bp.len());
            if max_len > 0 && dist <= max_edit_dist && dist * 100 / max_len <= 20 {
                canon_to_family.insert(other_canon.clone(), family_id);
            }
        }

        family_id += 1;
    }

    // Assign family IDs to arrays
    for (canon, indices) in &canon_groups {
        let fid = canon_to_family[canon];
        for &idx in indices {
            arrays[idx].family_id = Some(fid);
        }
    }

    eprintln!("Clustered into {} families", family_id);
}

/// Enrich: find smaller arrays (<min_length) that match known family consensuses
fn enrich_from_small(
    monomers_path: &str,
    families: &[(usize, Vec<u8>)],  // (family_id, canonical consensus)
    min_length: usize,
    max_edit_dist: usize,
    out: &mut BufWriter<File>,
) -> io::Result<usize> {
    let file = File::open(monomers_path)?;
    let reader = BufReader::new(file);
    let mut header_map: HashMap<String, usize> = HashMap::new();
    let mut enriched = 0usize;
    let mut seen: HashMap<String, bool> = HashMap::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        if line_num == 0 {
            for (i, field) in line.split('\t').enumerate() {
                header_map.insert(field.to_string(), i);
            }
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        let get = |name: &str| -> &str {
            header_map.get(name).and_then(|&i| fields.get(i)).copied().unwrap_or("")
        };

        let aid = get("array_id").to_string();
        if seen.contains_key(&aid) { continue; }

        let rtype = get("type");
        if rtype != "monomer" { continue; }

        let parts: Vec<&str> = aid.split('_').collect();
        if parts.len() < 5 { continue; }
        let arr_len: usize = parts[3].parse().unwrap_or(0);
        // Only look at arrays smaller than min_length (the ones not already classified)
        if arr_len >= min_length || arr_len < 100 { continue; }

        let cut_seq = get("cut_sequence").as_bytes().to_vec();
        let upper: Vec<u8> = cut_seq.iter().map(|b| b.to_ascii_uppercase()).collect();
        if upper.is_empty() { continue; }
        let base = find_base_period(&upper);

        // Try to match against known families
        for (fid, fam_canon) in families {
            let len_ratio = base.len().max(fam_canon.len()) as f64 / base.len().min(fam_canon.len()).max(1) as f64;
            if len_ratio > 2.0 { continue; }

            let dist = min_rotation_distance(&base, fam_canon);
            let max_len = base.len().max(fam_canon.len());
            if max_len > 0 && dist <= max_edit_dist && dist * 100 / max_len <= 20 {
                writeln!(out, "{}\t{}\t{}\t{}\t{}\tSF{:04}\tenriched",
                    parts[0], parts[1], parts[2], arr_len, parts[4], fid)?;
                enriched += 1;
                seen.insert(aid.clone(), true);
                break;
            }
        }
        seen.insert(aid, true);
    }

    Ok(enriched)
}

// ── Main ───────────────────────────────────────────────────────────────

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 || args.len() > 5 {
        eprintln!("Usage: sat-family <monomers.tsv> <output.tsv> [min_array_kb] [max_edit_dist]");
        eprintln!("");
        eprintln!("  Clusters satellite DNA arrays into families using monomer consensus");
        eprintln!("  with cyclic rotation, strand, and periodicity normalization.");
        eprintln!("");
        eprintln!("  Step 1: Cluster arrays >= min_array_kb into families");
        eprintln!("  Step 2: Enrich with smaller fragments matching family consensuses");
        eprintln!("");
        eprintln!("  min_array_kb:   minimum array length for initial clustering (default: 10)");
        eprintln!("  max_edit_dist:  maximum edit distance for merging families (default: 3)");
        process::exit(1);
    }

    let monomers_path = &args[1];
    let output_path = &args[2];
    let min_kb: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(10);
    let max_ed: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(3);
    let min_length = min_kb * 1000;

    // Step 1: Load large arrays
    eprintln!("Loading arrays >= {}kb...", min_kb);
    let mut arrays = load_arrays(monomers_path, min_length).unwrap_or_else(|e| {
        eprintln!("Error loading {}: {}", monomers_path, e);
        process::exit(1);
    });
    eprintln!("Loaded {} arrays", arrays.len());

    // Step 2: Cluster into families
    eprintln!("Clustering by canonical monomer consensus...");
    cluster_families(&mut arrays, max_ed);

    // Collect family consensuses for enrichment
    let mut family_consensus: HashMap<usize, Vec<u8>> = HashMap::new();
    for arr in &arrays {
        if let Some(fid) = arr.family_id {
            family_consensus.entry(fid).or_insert_with(|| arr.canonical.clone());
        }
    }

    // Step 3: Write output
    let out_file = File::create(output_path).unwrap_or_else(|e| {
        eprintln!("Error creating {}: {}", output_path, e);
        process::exit(1);
    });
    let mut out = BufWriter::new(out_file);

    writeln!(out, "# sat-family v0.1.0").unwrap();
    writeln!(out, "# min_array: {}kb  max_edit_dist: {}", min_kb, max_ed).unwrap();
    writeln!(out, "chr\tstart\tend\tarray_length\tperiod\tfamily\tconsensus\tbase_period\tsource").unwrap();

    let mut family_counts: HashMap<usize, usize> = HashMap::new();
    for arr in &arrays {
        let fid = arr.family_id.unwrap_or(9999);
        *family_counts.entry(fid).or_insert(0) += 1;
        let consensus_str = String::from_utf8_lossy(&arr.consensus);
        let base_str = String::from_utf8_lossy(&arr.base_period);
        writeln!(out, "{}\t{}\t{}\t{}\t{}\tSF{:04}\t{}\t{}\tcore",
            arr.chr, arr.start, arr.end, arr.end - arr.start,
            arr.period, fid, consensus_str, base_str).unwrap();
    }

    // Step 4: Enrich from small arrays
    eprintln!("Enriching from smaller arrays...");
    let fam_list: Vec<(usize, Vec<u8>)> = family_consensus.into_iter().collect();
    let enriched = enrich_from_small(monomers_path, &fam_list, min_length, max_ed, &mut out)
        .unwrap_or(0);

    out.flush().unwrap();

    // Summary
    eprintln!("\nResults:");
    eprintln!("  Core arrays (>= {}kb): {}", min_kb, arrays.len());
    eprintln!("  Families: {}", family_counts.len());
    eprintln!("  Enriched fragments: {}", enriched);

    // Top families
    let mut sorted_fams: Vec<_> = family_counts.iter().collect();
    sorted_fams.sort_by(|a, b| b.1.cmp(a.1));
    eprintln!("\n  Top families:");
    for (fid, count) in sorted_fams.iter().take(15) {
        let arr = arrays.iter().find(|a| a.family_id == Some(**fid)).unwrap();
        let base_str = String::from_utf8_lossy(&arr.base_period);
        eprintln!("    SF{:04}: {:>4} arrays, base_period={}bp ({})",
            fid, count, arr.base_period.len(), base_str);
    }
}
