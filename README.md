# unique_kmer_counter
Count number of unique kmers from fasta or fasta.gz files

This extremely simple tool counts the exact number of unique k-mers from a (multi)-fasta or a (multi)-fasta.gz file contains. 

- Only kmers on the A,C,G,T alphabet are counted. Lowercase letters (a,c,g,t) are considered as (A,C,G,T).
- Only kmers of size <= 32 are counted
- No canonicalisation
- No differentiation between sequences. If the input file contains more than a sequence (reads, chromosomes) they are all considered together, but not concatenated (no creation of alien kmers)

# Install
- clone: `git clone https://github.com/pierrepeterlongo/unique_kmer_counter`
- compile: `cd unique_kmer_counter && RUSTFLAGS="-C target-cpu=native" cargo install --path .`

# Usage 
`unique_kmer_counter <k> <input_file>`

# (big) Example
- Get the hg38 human genome: 
  - `wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz`
- Run (requires 32G of RAM):
  - `unique_kmer_counter 27 hg38.fa.gz`
- The program outputs intermediate results and ends with 
  ```
  Total nucleotides: 3209286105
  Number of distinct 27-mers: 2490500607
  ```

This requires 32G or RAM. 

# TODO and LIMITATIONS
The program was written in a few minutes. But, as I did not find any equivalent, I'm happy to share it here. 
However, I coded it for kmers of length <=32 (coded on 64 bits each) and I reserved 3 billions entries in a hash set as I used it for human genomes. 
- [X] Check options & use clap
- [ ] Adapt coding to kmer size
- [ ] Use also fastq[.gz] as input
- [ ] Print more stats



