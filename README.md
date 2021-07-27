# KARGA_Housekeeping
Multi-platform Toolkit for Antibiotic Resistance Identification in Housekeeping Genes from High-Throughput Sequencing Data.


# Installation
KARGA_Housekeeping requires the Java Virtual Machine (https://www.java.com/en/). The .class files available on this GitHub have been compiled on MS Windows 10 using 64-bit javac v.15.

# Usage
- KARGA_Housekeeping is based on k-mer probabilistic matching, similarly to his sibling KARGA (https://github.com/DataIntellSystLab/KARGA)
- KARGA_HouseKeeping can be launched from the command line. The minimum input is a read file in (optionally gzipped) FASTQ format, which is automatically detected if the extension is .fastq or .gz. 
- Without other parameters, our manually curated housekeeping genes database (housekeeping_db.fasta) is used with a default value of k=6. Please download the latest MEGARes release here: (curl) https://megares.meglab.org/download/index.php; https://megares.meglab.org/download/megares_v2.00/megares_full_database_v2.00.fasta.
- If you wish to classify mobile genetic elements (MGEs), we recommend to enable the multinomial classification option "m:y" which outputs multiple weighted hits. Several MGE reference databases can be used, e.g. ICEBerg (https://db-mml.sjtu.edu.cn/ICEberg/; curl https://db-mml.sjtu.edu.cn/ICEberg2/download/ICE_seq_all.fas).
- By default, the program outputs individual read classification as well as mapping of the resistome database given in input.
- If you want to classify antibiotic resistance in housekeeping genes, we recommend to use the dedicated KARGA_Housekeeping module (link TBA).
- Type "java KARGA readfile.fastq" for the default execution.
The java class accepts the following optional parameters: "k:your_k_value" (positive integer for k-mer length); "d:your_db_fasta" (any ARG/MGE database in FASTA format where resistance annotation is specified in the header); "f:your_read_fastq" (read file in FASTQ format with any file extension); "r:[y,yes,n,no]" (if you want to print or omit individual read classification, as the program is slightly faster when this print is omitted); "m:[y,yes,n,no]" (if you want to print the top-scoring ARG hits for each read, not only the best one); "i:your_value" (number of iterations to calculate frequency threshold from customized random string hit distribution, default is 125,000). For large databases, we recommend to use -Xmx16GB or larger.

# Output
- inputFileName_KARGA_mappedGenes.csv : a CSV file --one line per read-- with the following fields: Read_Idx, GeneProbability/KmersHitsOnGene/KmersHitsOnAllGenes/KmersTotal, GeneAnnotation.
- inputFileName_KARGA_mappedReads.csv : a CSV --one line per ARG-- with the following fields: GeneIdx, PercentGeneCovered, AverageKMerDepth. Note that ARGs with coverage below 1% are not printed.

