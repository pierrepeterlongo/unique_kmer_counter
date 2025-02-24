use clap::{Arg, Command};
use dashmap::DashSet;
use flate2::read::GzDecoder;
use rayon::ThreadPoolBuilder;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::process;
use std::sync::Arc;

const fn nucleotide_to_bits(n: u8) -> Option<u64> {
    match n {
        b'A' => Some(0b00),
        b'C' => Some(0b01),
        b'G' => Some(0b10),
        b'T' => Some(0b11),
        _ => None,
    }
}

fn kmer_to_u64(sequence: &[u8]) -> Option<u64> {
    if sequence.len() > 32 {
        return None;
    }

    let mut encoded: u64 = 0;
    for &nucleotide in sequence {
        encoded = (encoded << 2) | nucleotide_to_bits(nucleotide)?;
    }
    Some(encoded)
}

fn get_reader(filename: &str) -> io::Result<Box<dyn BufRead>> {
    let path = Path::new(filename);
    let file = File::open(path)?;

    if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn process_fasta_parallel(filename: &str, k: usize, reserve_size: usize) -> io::Result<(usize, usize, usize, usize)> {
    let reader = get_reader(filename)?;
    let lines: Vec<String> = reader.lines().filter_map(Result::ok).collect();

    let kmers = Arc::new(DashSet::with_capacity(reserve_size));
    let total_nucleotides = Arc::new(AtomicUsize::new(0));
    let nb_total_kmers = Arc::new(AtomicUsize::new(0));
    let nb_valid_kmers = Arc::new(AtomicUsize::new(0));

    lines.rchunks(1000).for_each(|chunk| {
        let local_kmers = DashSet::new();
        let mut sequence = Vec::new();
        let mut local_nuc_count = 0;
        let mut local_total_kmers = 0;
        let mut local_valid_kmers = 0;

        for line in chunk {
            if line.starts_with('>') {
                if !sequence.is_empty() {
                    local_total_kmers += sequence.len().saturating_sub(k) + 1;
                    for window in sequence.windows(k) {
                        if !window.contains(&b'N') {
                            if let Some(compact_kmer) = kmer_to_u64(window) {
                                local_kmers.insert(compact_kmer);
                                local_valid_kmers += 1;
                            }
                        }
                    }
                    sequence.clear();
                }
            } else {
                let clean_line = line.trim().bytes().map(|b| b.to_ascii_uppercase());
                local_nuc_count += line.trim().len();
                sequence.extend(clean_line);
            }
        }

        // Process last sequence
        if !sequence.is_empty() {
            local_total_kmers += sequence.len().saturating_sub(k) + 1;
            for window in sequence.windows(k) {
                if !window.contains(&b'N') {
                    if let Some(compact_kmer) = kmer_to_u64(window) {
                        local_kmers.insert(compact_kmer);
                        local_valid_kmers += 1;
                    }
                }
            }
        }

        // Merge results into shared atomic values
        total_nucleotides.fetch_add(local_nuc_count, Ordering::Relaxed);
        nb_total_kmers.fetch_add(local_total_kmers, Ordering::Relaxed);
        nb_valid_kmers.fetch_add(local_valid_kmers, Ordering::Relaxed);

        for kmer in local_kmers.iter() {
            kmers.insert(*kmer);
        }
    });

    Ok((
        kmers.len(),
        total_nucleotides.load(Ordering::Relaxed),
        nb_total_kmers.load(Ordering::Relaxed),
        nb_valid_kmers.load(Ordering::Relaxed),
    ))
}



fn process_fasta_parallel_only_count(filename: &str, k: usize) -> io::Result<(usize, usize, usize)> {
    let reader = get_reader(filename)?;
    let lines: Vec<String> = reader.lines().filter_map(Result::ok).collect();

    let total_nucleotides = Arc::new(AtomicUsize::new(0));
    let nb_total_kmers = Arc::new(AtomicUsize::new(0));
    let nb_valid_kmers = Arc::new(AtomicUsize::new(0));

    lines.rchunks(1000).for_each(|chunk| {
        let mut sequence = Vec::new();
        let mut local_nuc_count = 0;
        let mut local_total_kmers = 0;
        let mut local_valid_kmers = 0;

        for line in chunk {
            if line.starts_with('>') {
                if !sequence.is_empty() {
                    local_total_kmers += sequence.len().saturating_sub(k) + 1;
                    for window in sequence.windows(k) {
                        if !window.contains(&b'N') {
                            if let Some(_) = kmer_to_u64(window) {
                                local_valid_kmers += 1;
                            }
                        }
                    }
                    sequence.clear();
                }
            } else {
                let clean_line = line.trim().bytes().map(|b| b.to_ascii_uppercase());
                local_nuc_count += line.trim().len();
                sequence.extend(clean_line);
            }
        }

        // Process last sequence
        if !sequence.is_empty() {
            local_total_kmers += sequence.len().saturating_sub(k) + 1;
            for window in sequence.windows(k) {
                if !window.contains(&b'N') {
                    if let Some(_) = kmer_to_u64(window) {
                        local_valid_kmers += 1;
                    }
                }
            }
        }

        // Merge results into shared atomic values
        total_nucleotides.fetch_add(local_nuc_count, Ordering::Relaxed);
        nb_total_kmers.fetch_add(local_total_kmers, Ordering::Relaxed);
        nb_valid_kmers.fetch_add(local_valid_kmers, Ordering::Relaxed);
    });

    Ok((
        total_nucleotides.load(Ordering::Relaxed),
        nb_total_kmers.load(Ordering::Relaxed),
        nb_valid_kmers.load(Ordering::Relaxed),
    ))
}


fn main() {
    let matches = Command::new("Unique Kmer Counter")
        .version("1.0")
        .author("Author Name <author@example.com>")
        .about("Counts unique kmers in a FASTA file")
        .arg(
            Arg::new("k")
                .short('k')
                .long("kmer-size")
                .value_name("K")
                .help("Sets the k-mer size")
                .required(true)
                .num_args(1),
        )
        .arg(
            Arg::new("fasta_file")
                .short('f')
                .long("input-file")
                .help("Sets the input FASTA file")
                .required(true)
                .num_args(1),
        )
        .arg(
            Arg::new("reserve_size")
                .short('r')
                .long("reserve")
                .value_name("RESERVE")
                .help("Sets the initial reserve size for the HashSet. \
                Useless with the only_count option")
                .default_value("3000000000")
                .num_args(1),
        )
        .arg(
            Arg::new("only_count")
            .short('c')
            .long("only-count")
            .num_args(0) 
            .help("Only count the number of kmers and nucleotides")
        )
        .arg(
            Arg::new("max_threads")
                .short('t')
                .long("max-threads")
                .value_name("THREADS")
                .help("Limits the maximum number of threads")
                .default_value("0")
                .num_args(1),
        )
        .get_matches();

        let k = matches
            .get_one::<String>("k")
            .and_then(|s| s.parse::<usize>().ok())  // Parse safely
            .unwrap_or_else(|| {
                eprintln!("Error: k must be a positive integer");
                process::exit(1);
        });


    if k > 32 {
        eprintln!("Error: k must be less than or equal to 32");
        process::exit(1);
    }

    if k < 1 {
        eprintln!("Error: k must be positive");
        process::exit(1);
    }

    let fasta_file = matches.get_one::<String>("fasta_file").unwrap();
    let reserve_size = matches
        .get_one::<String>("reserve_size")
        .and_then(|s| s.parse::<usize>().ok())  // Parse safely
        .unwrap_or_else(|| {
            eprintln!("Error: reserve_size must be a positive integer");
            process::exit(1);
        });
    
    let max_threads = matches
        .get_one::<String>("max_threads")
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0);

    if max_threads > 0 {
        ThreadPoolBuilder::new()
            .num_threads(max_threads)
            .build_global()
            .expect("Failed to set thread limit");
    }
    if matches.get_flag("only_count") {
        match process_fasta_parallel_only_count(fasta_file, k) {
            Ok((nuc_count, nb_total_kmers, nb_valid_kmers)) => {
                println!("Total nucleotides: {}", nuc_count);
                println!("Total k-mers: {}", nb_total_kmers);
                println!("Valid k-mers: {}", nb_valid_kmers);
            }
            Err(e) => {
                eprintln!("Error processing file: {}", e);
                process::exit(1);
            }
        }
        return;
    }

    match process_fasta_parallel(fasta_file, k, reserve_size) {
        Ok((kmer_count, nuc_count, nb_total_kmers, nb_valid_kmers)) => {
            println!("Total nucleotides: {}", nuc_count);
            println!("Total k-mers: {}", nb_total_kmers);
            println!("Valid k-mers: {}", nb_valid_kmers);
            println!("Number of distinct {}-mers: {}", k, kmer_count);
        }
        Err(e) => {
            eprintln!("Error processing file: {}", e);
            process::exit(1);
        }
    }
}

