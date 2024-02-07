# Kmer-Counter is a tool for gathering data about bacteria

This is a tool to count the amount of times a k-mer appears across multiple fna/fasta files. It also keeps track of how many of the seen kmers are resistant or susceptible(to an antibiotic).

## Example

Say we have two files A and B. A is classified as resistant and B as susceptible. Suppose we come across these 4-mers. ACCA, GTCA, GGAC. If the first 4-mer is present only in A, the second 4-mer in only B and GGAC in both the program would write down the following:

4-mer, res, sus
20, 1, 0
180, 0, 1
161, 1, 1

Note that the only thing we're counting is the kmers presence not how many times it appears in the file.

The kmers are encoded as numbers to make calculations faster (a=00, c=01, g=10, t=11)

## Usage

**UPON RUNNING THE SCRIPT FOR THE FIRST TIME, THE USER IS PROMPTED FOR THE ABSOLUTE PATH OF THE FOLDER WHERE THE FASTA FILES ARE STORED. WE RECOMMEND MAKING A FOLDER THAT HOUSES FOLDERS THAT STORE THE FASTA FILES.**

There are two ways the script can be run. Once the script finishes its output will be written to a file called "counts.csv"

### 1

```bash
./kmerCounter
```

The program asks the user to input the required parameters instead of giving them as command line arguments like in the next heading

### 2

```bash
./kmerCounter folder threads k
```

- 'folder' is the folder where the fna files are stored
- 'threads' is how many threads you want the program to use(Physical threads/cores)
- k is the length of the kmer(4-mer = aaaa, 2-mer = aa etc)
