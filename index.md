# Detailed Bioinformatics pipelines for Short- and Long-read sequences with Kraken2

*- Funded by [European COST Action](https://www.cost.eu) and [EuroMIC](https://www.euro-mic.org)*

---

## Summary

This guide outlines how Short- and Long-read sequences can be analysed using a custom pipeline
to assess data from 16S rRNA sequencing of microbial communities associated with MIC
(Microbiologically Influenced Corrosion). 

It is designed to be accessible to readers with limited
technical expertise, ensuring clarity and ease of understanding. Terminal commands and their
parameters will be explained objectively. When applicable, other, openly available tools that are suited
for the task, is listed as an alternative option.

### Pre-requisites

- Workstation with at least ***4 CPU Cores***, ***16GB RAM*** and ***75GB of free space*** (8+ Cores, 32GB+ RAM and 150GB free space recommended).
- Basic understanding of CLI (Command Line Interface) to use the listed tools
- *Optional*: Knowledge of Bash or other scripting/programming language (Python, R, etc.)


### Terminal commands

CLI or Terminal commands are in Bash (which is the default language of Linux terminals). Many of them uses placeholder reference to a data or location `<some/location>`. This placeholder will not work, as it is, on your computer, as long as does not refer to a real source.


*For example:*

```bash
kraken2-build --build --db <database/location>
```

- Make sure to replace the placeholder with direct path, otherwise the command will result in error. 
- Remove the `< >` symbols for the correct path.

Some commands are with specific parameters, for example: `-q 30`; setting depends on the project, dataset and requirements, make sure you do your research what is the optimal parameter for your specific use case. If you are not sure, you may use the listed ones, but keep in mind that may result incorrect or misleading results.


### Multiple files handling

Most likely you will have a dataset with many-many files. As the commands are written for single use, optimally you want to loop through all the files in order to generate the result from all of them, also making sure the output file naming does not overwrite any existing and important file.

While the guide should work without any programming knowledge, for time-saving I highly recommend to learn the concept of Loops. Which is basically just jumping to each file and apply a command them one-by-one quickly.

*Bash loop example:*

```text
Dataset (paired-end):

SEQ_01_R1.fastq
SEQ_01_R2.fastq
SEQ_02_R1.fastq
SEQ_02_R2.fastq
SEQ_03_R1.fastq
SEQ_03_R2.fastq
SEQ_04_R1.fastq
SEQ_04_R2.fastq
SEQ_05_R1.fastq
SEQ_05_R2.fastq
```

```bash
#!/bin/bash

# Loop through files from SEQ_01 to SEQ_05
for i in {01..05}; do
    R1="SEQ_${i}_R1.fastq"
    R2="SEQ_${i}_R2.fastq"
    
    echo "Processing $R1 and $R2..."
    
    fastqc -o SEQ_${i} "$R1" "$R2"
done
```
The script will read the names of the files, running FastQC with pair-end input and creates a result with the pair name.

If the concept seems difficult here are some brief tutorials:
- [Bash in general](https://www.freecodecamp.org/news/bash-scripting-tutorial-linux-shell-script-and-command-line-for-beginners/)
- [Loops in Bash](https://www.geeksforgeeks.org/bash-scripting-for-loop/)

### Windows

If a Linux or Mac system is not available, the easiest way for Windows users is the [Git for Windows](https://gitforwindows.org) tool, that includes a Git Bash application that behaves almost identical to a native Unix command line (with some limitations).

Path: use `/c/User/YourFolder` to access `C:\User\YourFolder`

Running Bash file:

- Open Notepad and paste the following: 
```text
#!/bin/bash

echo "Hello from Git Bash"
```
- Save it as `my_script.sh` somewhere you will find. Navigate to that location, right click on an empty space, if Git for Windows is installed, you should see in the dropdown, open it and run: `bash my_script.sh` If you see the `Hello from Git Bash` message, it works as expected.

---

## Short reads - Data Preparation

### Goal:

Quality Control (QC) and trimming of short reads are essential to ensure data accuracy and reliability by addressing sequencing errors, adapter contamination, and low-quality bases, which can lead to incorrect taxonomic assignments or misinterpretation of microbial diversity. These steps improve downstream analyses, such as taxonomic classification. Additionally, QC and trimming reduce computational burden by decreasing sequence size, making downstream processes computationally less heavy.

### **1. Obtain Raw Sequencing Data**
- **Format:** Paired-end (`reads_R1.fastq` and `reads_R2.fastq`) or single-end (`reads.fastq`).
- **Source:** Raw data from sequencing platforms (Illumina, etc.).

### **2. Quality Control (Pre-Trim)**
- **Tool:** **[FastQC](https://github.com/s-andrews/FastQC)** to assess read quality.
  ```bash
  fastqc -o reads_fastqc <reads_R1.fastq> <reads_R2.fastq>
  ```
    **Parameter:**
    `-o`: Output file name

    FastQC generates several files that not necessarily useful independently. To summarize the output files, use **[MultiQC](https://github.com/MultiQC/MultiQC)** to merge the FastQC files and analyze.

    ```bash
    multiqc <FASTQC_DIR> -o <REPORT_DIR>
    ```
  Navigate to the folder that set as OUTPUT_DIR and find a `multiqc_report.html` for details and visualizations.

### **3. Trimming**
- **Tool:** **[Cutadapt](https://github.com/marcelm/cutadapt)**, Trimmomatic, BBDuk or equivalent to clean reads.

  ```bash
  cutadapt -a <FORWARD_ADAPTER_SEQ> -A <REVERSE_ADAPTER_SEQ> \
  -q 20 -O 5 -m 36 \
  -o trimmed_R1.fastq -p trimmed_R2.fastq \
  <reads_R1.fastq> <reads_R2.fastq>
  ```
  Common adapter sequences for Illumina 16S amplicon sequencing are:

  **Forward adapter:**
  
  ```
  AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
  ```

  **Reverse adapter:**
  ```
  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
  ```
  
    **Parameters:**

    `-a`: Forward adapter
 
    `-A`: Reverse adapter
  
    `-q`: Quality cutoff (Phred score threshold) - below are trimmed from the 3' and/or 5' ends of the reads
    
    `-O`: Minimum overlap length - minimum number of overlapping bases required between the read and the adapter for the adapter to be trimmed
    
    `-m`: Minimum length - minimum allowed length for a read after trimming
    
    `-o`: Output file name

### **4. Post-Process Quality Control**
- Run **FastQC** and **MultiQC** again on the trimmed files to ensure the quality is sufficient for downstream analysis.

    **You can compare the Raw and Trimmed quality**:
    ```bash
    multiqc <FASTQC_PRE_TRIM_DIR> <FASTQC_POST_TRIM_DIR> -o <REPORT_DIR>
    ```
---

## Long reads - Data Preparation

### Goal:

Quality control (QC), filtering, and trimming are essential for processing long reads, like Oxford Nanopore (ONT) sequences to ensure accuracy, efficiency, and usability in genomic analysis. Nanopore sequencing has a higher error rate than Illumina, including insertions, deletions, and mismatches, making QC and trimming crucial for reducing errors. Adapters used during library preparation may remain in the reads and need to be removed to avoid interference with downstream analyses. These steps further remove low-quality bases, and chimeric reads, improving genome assembly, variant calling, and functional annotation. Additionally, filtering out low-quality and excessively short reads reduces computational burden and enhances alignment accuracy.

### **1. Obtain Raw Sequencing Data**
- **Format:** Single-end (`reads.fastq`).
- **Source:** Raw data from sequencing platforms (Nanopore, etc.).

### **2. Quality Control, Trim and Filter**
   - Tool: **[FastPLong](https://github.com/OpenGene/fastplong)**
   - Purpose: Evaluate read length distribution, quality scores, and sequencing statistics.

        ```bash
        fastplong -i <reads.fastq> -o <OUTPUT_FILENAME> -h QUALITY_CONTROL.html -j QUALITY_CONTROL.json -l 320
        ```
        **Parameters:**

        `-i`: Input file
 
        `-o`: Output file name
  
        `-h`: Quality control report .html file
    
        `-j`: Quality control report .json file

        `-l`: Reads shorter than set number, will be discarded


- FastPLong generates quality control files as well (html and json). To summarize the output files, use **[MultiQC](https://github.com/MultiQC/MultiQC)** to summarizing the reports.

    ```bash
    multiqc <QC_DIR> -o <REPORT_DIR>
    ```
  Navigate to the folder that set as OUTPUT_DIR and find a `multiqc_report.html` for details and visualizations.

---
From this point the steps can be applied to Short- and Long-reads sequences too.

---

## Taxonomic Classification

### Goal:
Taxonomic classification of DNA or RNA sequences used to find exact matches in pre-built database of reference genomes.

Kraken2 using exact k-mer matches to compare sequences against a reference database, assigning them to specific taxonomic groups with high accuracy. Kraken2 is particularly useful for metagenomic analysis, enabling the identification and quantification of microbial communities in complex samples. It is widely used in clinical diagnostics for pathogen detection and in environmental studies for microbiome profiling. Its flexibility, speed, and ability to handle both short and long reads make it a valuable tool for researchers.


### **1. Install Kraken2**
- Follow the instructions to install Kraken2 from [Kraken2 GitHub](https://github.com/DerrickWood/kraken2)

- See [more details](https://ccb.jhu.edu/software/kraken2/) and publications about Kraken2.
- Kraken2 needs to be installed in its own repository folder, then scripts can be copied to another location if prefered.
- Scripts: `kraken2`, `kraken2-build`, `kraken2-inspect`

### **2. Build the Kraken2 Database**
This guide uses SILVA-138 SSURef-NR99 database, but Kraken2 has several other built-in databases and able to create from [external sources](https://benlangmead.github.io/aws-indexes/k2).

- **Create a database directory:**
   ```bash
   mkdir /kraken2_silva
   ```
- **Add SILVA sequences to the Kraken2 database:**
   ```bash
   kraken2-build --special silva --db /kraken2_silva
   ```

   - **Note:** This process may rbe **memory (RAM) intensive** and sufficient disk space necessary, depending on the database size.

### **3. Run Kraken2 for Taxonomic Classification**


- **Paired-end reads (Short):**
  ```bash
  kraken2 --db </path/to/kraken2_silva> \
    --paired <trimmed_R1.fastq> <trimmed_R2.fastq> \
    --threads <your_cpu> \
    --confidence 0.2 \
    --output kraken2_output.txt \
    --report kraken2_report.tabular
  ```
- **Single-end reads (Long):**
  ```bash
  kraken2 --db </path/to/kraken2_silva> \
    <trimmed_reads.fastq> \
    --threads <your_cpu> \
    --confidence 0.2 \
    --output kraken2_output.txt \
    --report kraken2_report.tabular
  ```
    
    **Parameters:**
    
  `--db`: Built Kraken2 database location
    
  `--paired`: Paired-end reads, usually short-read sequences are paired-end
 
  `--threads`: Your CPU threads, if not sure, use `8`
  
  `--confidence`: Confidence score threshold must be between 0 and 1, see more details below
    
  `--output`: Read sequence output

  `--report`: Taxonomic result, file name can be anything readable (.txt, etc.)

  
- **Confidence score**

    Confidence score is a numerical value that represents the reliability or certainty of a prediction, classification, or analysis result. It provides a measure of how likely a computational result is to be correct, helping researchers prioritize high-confidence findings and filter out uncertain ones.
  
    - Articles:
  
      - [Impact of database choice and confidence score on the performance of taxonomic classification using Kraken2](https://link.springer.com/article/10.1007/s42994-024-00178-0)
      - [An Alignment Confidence Score Capturing Robustness to Guide Tree Uncertainty](https://pmc.ncbi.nlm.nih.gov/articles/PMC2908709/)
      - [The standardisation of the approach to metagenomic human gut analysis: from sample collection to microbiome profiling](https://www.nature.com/articles/s41598-022-12037-3)

    
**IMPORTANT**: If you encounter a result like `100% classification`, it is a **false positive**, which is likely caused by incorrect database.





----
Seqkit
Bracken
Visualization