# KARGVA
K-mer-based Antibiotic Resistance Gene Variant Analyzer (KARGVA). KARGVA is a Multi-platform Toolkit for Identification of Antibiotic Resistance from Sequencing Data Conferred by Point Mutations in Bacterial Genes.

# Installation
KARGVA requires the Java Virtual Machine (https://www.java.com/en/). The .class files available on this GitHub have been compiled on MS Windows 10 using 64-bit javac v.15.

# Usage
- KARGVA is based on k-mer probabilistic matching, similarly to his sibling KARGA (https://github.com/DataIntellSystLab/KARGA)
- KARGVA can be launched from the command line. The minimum input is a read file in (optionally gzipped) FASTQ format, which is automatically detected if the extension is .fastq or .gz. 
- Without other parameters, our manually curated database is used (the latest update is kargva_db_v5.fasta). Our database includes antibiotic resistance gene variants (ARGVs), where the antibiotic resistance is due to specific point mutations in a gene (e.g. in housekeeping genes), integrated from different sources, namely: MEGARes (https://megares.meglab.org/); CARD (https://card.mcmaster.ca/); NDARO (https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/). We eliminated sequence redundancy, but kept some heterogeneity in the combination of point mutations. Thus, each sequence can be traced back to one or more original sources; the header contains also information about the combination of aminoacid mutations conferring antibiotic resistance included in the given sequence.
- Type "java KARGVA readfile.fastq" for the default execution.
The java class accepts the following optional parameters: "k:your_k_value" (positive integer for k-mer length, as aminoacids, default value is 9); "d:your_db_fasta" (any database of housekeeping genes conferring antibiotic resistance in FASTA format, where the resistance mutations are specified in the header); "f:your_reads_fastq" (read file in FASTQ format with any file extension); "i:your_value" (number of iterations to calculate frequency threshold from customized random string hit distribution, default is 12,500); "m:[y,yes,n,no]" (if you want to print the top-scoring ARGV hits for each read, not only the best one). For large databases, we recommend to use -Xmx16GB or larger.

# Output
- inputFileName_KARGVA_mappedReads.csv : a CSV file --one line per read-- with the following fields: Read_Idx, GeneAnnotation, GeneScore/KmerSNPsHits/KmersHitsOnGene/KmersHitsOnAllGenes/KmersTotal. Given that antibiotic resistance in housekeeping genes can be due to one or more aminoacidic mutations, in different combinations, we report all the best scoring genes and point mutations present in the database within 5% tolerance ratio compared to the top-scoring hit.
- inputFileName_KARGVA_mappedGenes.csv : a CSV --one line per ARG-- with the following fields: GeneIdx, PercentGeneCovered, AverageKMerDepth. Note that genes with coverage below 1% are not printed.

# References
Preprint paper: https://www.biorxiv.org/content/10.1101/2022.08.12.503773v1
