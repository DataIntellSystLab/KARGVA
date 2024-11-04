# KARGVA
K-mer-based Antibiotic Resistance Gene Variant Analyzer (KARGVA). KARGVA is a Multi-platform Toolkit for Identification of Antibiotic Resistance from Sequencing Data Conferred by Point Mutations in Bacterial Genes.

# Installation
KARGVA requires the Java Virtual Machine (https://www.java.com/en/). The .class files available on this GitHub have been compiled on MS Windows 10 using 64-bit javac v.15.

# Usage
- KARGVA is based on k-mer probabilistic matching, similarly to his sibling KARGA (https://github.com/DataIntellSystLab/KARGA)
- KARGVA can be launched from the command line. The minimum input is a read file in (optionally gzipped) FASTQ format, which is automatically detected if the extension is .fastq or .gz. 
- Without other parameters, our manually curated database is used (the latest update is kargva_db_v5.fasta). Our database includes antibiotic resistance gene variants (ARGVs), where the antibiotic resistance is due to specific point mutations in a gene (e.g. in housekeeping genes), integrated from different sources, namely: MEGARes (https://megares.meglab.org/); CARD (https://card.mcmaster.ca/); NDARO (https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/). We eliminated sequence redundancy, but kept some heterogeneity in the combination of point mutations. Thus, each sequence can be traced back to one or more original sources; the header contains also information about the combination of aminoacid mutations conferring antibiotic resistance included in the given sequence.
- Type "java KARGVA readfile.fastq" for the default execution. For databases larger than the KARGVA default (kargva_db_v5.fasta), we recommend to use -Xmx16GB or larger.
- The java class accepts the following optional command line parameters: "k:your_k_value" (positive integer for k-mer length, as aminoacids, default value is 9, smaller values are better for higher error rates of the sequencer); "d:your_db_fasta" (any database of genes conferring antibiotic resistance thorugh mutations, in FASTA format, where the resistance mutations are specified in the header, see the kargva_db_v5.fasta example); "f:your_reads_fastq" (read file in FASTQ format with any file extension); "i:your_iterr_value (number of iterations to calculate false positive frequency threshold from a customized count hit distribution of random reads, default is 12,500, it will be calibrated on average read length and k-mer length, increase it if the data base size increases or the read length increases, but very large values can slow the program); "m:[y,yes,n,no]" (if you want to print the top-scoring ARGV hits for each read, not only the best one, default is yes). 

# Output
- inputFileName_KARGVA_mappedReads.csv : a CSV file --one line per read-- with the following fields: Read_Idx, GeneAnnotation, GeneScore/KmerSNPsHits/KmersHitsOnGene/KmersHitsOnAllGenes/KmersTotal. Given that antibiotic resistance in housekeeping genes can be due to one or more aminoacidic mutations, in different combinations, we report all the best scoring genes and point mutations present in the database within 5% tolerance ratio compared to the top-scoring hit.
- inputFileName_KARGVA_mappedGenes.csv : a CSV --one line per ARG-- with the following fields: GeneIdx, PercentGeneCovered, AverageKMerDepth. Note that genes with coverage below 1% are not printed.
- Regarding the output file fields, the mappedGenes.csv file include: GeneIdx - identifier of individual ARGV; PercentGeneCovered - percent of k-mers over the total number of distinct k-mers in the gene (weighting multiplicity); AverageKMerDepth- average number of times a k-mer is covered. The mappedReads.csv file include: Idx - identifier of individual read; GeneAnnotation - identifier of individual ARGV that the read maps to; GeneScore - score of ARGV mapping reliability; KmerSNPsHits - number of k-mers with resistance SNP found in the mapped ARGV; KmersHitsOnGene - number of k-mers found in the mapped ARGV; KmersHitsOnAllGenes - number of k-mers found in all ARGVs; KmersTotal - total number of k-mers in the read

# References
Preprint paper: https://www.biorxiv.org/content/10.1101/2022.08.12.503773v1
Published paper: https://pubmed.ncbi.nlm.nih.gov/36960290/

# Tutorial
Let's start by running KARGVA with options on our example file. To do so, open a terminal, go the directory where you cloned this repository, and run the following:

```java KARGVA k:17 testdata_simul.fastq```

(If you get an error based on the java runtime version, you can recompile KARGVA by running `javac KARGVA.java`.)

This command will run KARGVA on our default fastq file with k-mer length 17. If you have multiple fasta files that should be processed together, such as paired reads R1 and R2 files, you can concatenate them in a single fastq file.

## Mapped genes and reads
The first output file, `testdata_simul_KARGVA_mappedGenes.csv`, is the per-gene report. Each line of the file indicates a match with a gene. You can check the first lines by running the following: `head testdata_simul_KARGVA_mappedGenes.csv`.
- `GeneIdx`, is the identifier of individual ARGV, according to the database (header of the fasta match). With the KARGVA database, the `GeneIdx` is splitof into fields separated by `|`. For example: `>1002|V213I|ARO:3003901|Staphylococcus aureus GlpT with mutation conferring ...`. The first field is the KARGVA ARGV index, the second describes the variant at the protein level (wild type 'V' is changed into 'I'), and the others are the descriptors of CARD, MEGARes, and NDARO. Note that "NA" means that there is no corresponding entry for that field in the original database.
- `KmerSNPHits` is the number of k-mers with resistance SNP found in the mapped ARGV
- `PercentGeneCovered` is the percent of k-mers over the total number of distinct k-mers in the gene (weighting multiplicity)
- `AverageKMerDepth` is the average number of times a k-mer is covered

The second output file, `testdata_simul_KARGVA_mappedReads.csv`, is the per-read report. Each line of the file indicates a match with a read. You can check the first lines by running the following: `head testdata_simul_KARGVA_mappedReads.csv`.
- `Idx` is the read name, according to the fastq
- `GeneAnnotation` is the description (fasta header) of gene matched by the read. If you are using the KARGVA database, see above for a description of the fields.
- `GeneScore` is score of ARGV mapping reliability; `KmerSNPsHits` is the number of k-mers with resistance SNP found in the mapped ARGV; `KmersHitsOnGene` is the number of k-mers found in the mapped ARGV; `KmersHitsOnAllGenes` is the number of k-mers found in all ARGVs; and `KmersTotal` is total number of k-mers in the read.


