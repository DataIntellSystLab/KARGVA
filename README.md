# KARGA_Housekeeping
Multi-platform Toolkit for Antibiotic Resistance Identification in Housekeeping Genes from High-Throughput Sequencing Data.


# Installation
KARGA_Housekeeping requires the Java Virtual Machine (https://www.java.com/en/). The .class files available on this GitHub have been compiled on MS Windows 10 using 64-bit javac v.15.

# Usage
- KARGA_Housekeeping is based on k-mer probabilistic matching, similarly to his sibling KARGA (https://github.com/DataIntellSystLab/KARGA)
- KARGA_HouseKeeping can be launched from the command line. The minimum input is a read file in (optionally gzipped) FASTQ format, which is automatically detected if the extension is .fastq or .gz. 
- Without other parameters, our manually curated housekeeping genes database (housekeeping_db.fasta) is used. Our database includes bacterial housekeeping genes integrated from different sources, namely: MEGARes (https://megares.meglab.org/); CARD (https://card.mcmaster.ca/); NDARO (https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/). Each sequence can be traced back to one or more sources. The header contains also information about the aminoacid mutations conferring antibiotic resistance.
- Type "java KARGA_Housekeeping readfile.fastq" for the default execution.
The java class accepts the following optional parameters: "k:your_k_value" (positive integer for k-mer length, as aminoacids); "d:your_db_fasta" (any database of housekeeping genes conferring antibiotic resistance in FASTA format, where the resistance mutations are specified in the header); "f:your_read_fastq" (read file in FASTQ format with any file extension); "i:your_value" (number of iterations to calculate frequency threshold from customized random string hit distribution, default is 10,000); "m:[y,yes,n,no]" (if you want to print the top-scoring ARG hits for each read, not only the best one). For large databases, we recommend to use -Xmx16GB or larger.

# Output
- inputFileName_KARGAHK_mappedGenes.csv : a CSV file --one line per read-- with the following fields: Read_Idx, GeneAnnotation, GeneScore/KmerSNPsHits/KmersHitsOnGene/KmersHitsOnAllGenes/KmersTotal. Given that antibiotic resistance in housekeeping genes can be due to one or more aminoacidic mutations, in different combinations, we report all the best scoring genes and mutations present in the database within 10% tolerance compared to the top hit.
- inputFileName_KARGAHK_mappedReads.csv : a CSV --one line per ARG-- with the following fields: GeneIdx, PercentGeneCovered, AverageKMerDepth. Note that genes with coverage below 1% are not printed.

