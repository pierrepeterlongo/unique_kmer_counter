use clap::{Command, Arg};
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::process;
use flate2::read::GzDecoder;

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

fn count_kmers_and_nucleotides(filename: &str, k: usize, reserve_size: usize) -> io::Result<(usize, usize, usize, usize)> {
    let reader = get_reader(filename)?;
    let mut kmers = HashSet::new();
    kmers.reserve(reserve_size); 
    let mut sequence = Vec::new();
    let mut total_nucleotides = 0;
    let mut progress_counter = 0;
    let mut nb_total_kmers = 0;
    let mut nb_valid_kmers = 0;

    let progress_interval = 10_000_000; // Report every 1M nucleotides

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if sequence.len() !=0 {
                nb_total_kmers += sequence.len() - k + 1;
            }
            
            // Process previous sequence if any
            for window in sequence.windows(k) {
                if !window.contains(&b'N') {
                    if let Some(compact_kmer) = kmer_to_u64(window) {
                        kmers.insert(compact_kmer);
                        nb_valid_kmers += 1;
                    }
                }
            }
            sequence.clear();
        } else {
            let clean_line = line.trim().bytes().map(|b| b.to_ascii_uppercase());
            total_nucleotides += line.trim().len();
            progress_counter += line.trim().len();
            sequence.extend(clean_line);


            if progress_counter >= progress_interval {
                eprintln!("Processed {} million nucleotides, found {} kmers", 
                    total_nucleotides / 1_000_000, 
                    kmers.len());
                progress_counter = 0;
            }
        }
    }
    
    // Process the last sequence
    if sequence.len() !=0 {
        nb_total_kmers += sequence.len() - k + 1;
    }
    for window in sequence.windows(k) {
        if !window.contains(&b'N') {
            if let Some(compact_kmer) = kmer_to_u64(window) {
                kmers.insert(compact_kmer);
                nb_valid_kmers += 1;
            }
        }
    }

    Ok((kmers.len(), total_nucleotides, nb_total_kmers, nb_valid_kmers))
}

fn main() {
    let matches = Command::new("Unique Kmer Counter")
        .version("1.0")
        .author("Author Name <author@example.com>")
        .about("Counts unique kmers in a FASTA file")
        .arg(Arg::new("k")
            .short('k')
            .long("kmer-size")
            .value_name("K")
            .help("Sets the k-mer size")
            .required(true)
            .num_args(1))
        .arg(Arg::new("fasta_file")
            .short('f')
            .long("input-file")
            .help("Sets the input FASTA file")
            .required(true)
            .num_args(1))
        .arg(Arg::new("reserve_size")
            .short('r')
            .long("reserve")
            .value_name("RESERVE")
            .help("Sets the initial reserve size for the HashSet")
            .default_value("3000000000")
            .num_args(1))
        .get_matches();

    let k = matches.get_one::<String>("k").unwrap().parse::<usize>().unwrap_or_else(|_| {
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
    let reserve_size = matches.get_one::<String>("reserve_size").unwrap().parse::<usize>().unwrap_or_else(|_| {
        eprintln!("Error: reserve_size must be a positive integer");
        process::exit(1);
    });

    match count_kmers_and_nucleotides(fasta_file, k, reserve_size) {
        Ok((kmer_count, nuc_count, nb_total_kmers, nb_valid_kmers)) => {
            println!("Total nucleotides: {}", nuc_count);
            println!("Total k-mers: {}", nb_total_kmers);
            println!("Valid k-mers: {}", nb_valid_kmers);
            println!("Number of distinct {}-mers: {}", k, kmer_count);
        },
        Err(e) => {
            eprintln!("Error processing file: {}", e);
            process::exit(1);
        }
    }
}