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

*For example*:

```bash
kraken2-build --build --db <database/location>
```

- Make sure to replace the placeholder with direct path, otherwise the command will result in error. 
- Remove the `< >` symbols for the correct path.

Some commands are with specific parameters, for example: `-q 30`; setting depends on the project, dataset and requirements, make sure you do your research what is the optimal parameter for your specific use case. If you are not sure, you may use the listed ones, but keep in mind that may result incorrect or misleading results.

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
    **Parameter**:
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

  **Forward adapter**:
  
  ```
  AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
  ```

  **Reverse adapter**:
  ```
  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
  ```
  
    **Parameters**:

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
        **Parameters**:

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

----
(Both)
Kraken2 DB build
Taxonomic classification
Seqkit
Bracken
Visualization