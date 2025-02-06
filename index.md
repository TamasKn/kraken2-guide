# Detailed Bioinformatics pipelines for Short- and Long-read sequences with Kraken2

*- Author: Tamas Knisz* \
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
- *Optional*: Knowledge of Bash or other scripting/programming language (Python or R)


### Terminal commands

CLI or Terminal commands are in Bash (which is the default language of Linux terminals). Many of them use placeholder reference to a data or location `<some/location>`. This placeholder will not work, as it is, on your computer, as long as does not refer to an actual source.


*For example:*

```bash
kraken2-build --build --db <database/location>
```

- Make sure to replace the placeholder with direct path, otherwise the command will result in error. 
- **Remove** the `< >` symbols for the correct path.

Some commands are with specific parameters, for example: `-q 30`; setting depends on the project, dataset and requirements, make sure you do your research what is the optimal parameter for your specific use case. If you are not sure, you may use the listed or software default ones, but keep in mind that may lead to incorrect or misleading results.


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

If the concept seems difficult, here are some brief tutorials:
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
- Save it as `my_script.sh` somewhere and navigate to that location
- Right click on an empty space, if Git for Windows is installed, you should see in the dropdown, open it and run: `bash my_script.sh`
- If you see the `Hello from Git Bash` message, your terminal is ready to use

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

  FastQC generates several files that not necessarily useful independently. 


- To summarize the output files, use **[MultiQC](https://github.com/MultiQC/MultiQC)** to merge the FastQC files and analyze.

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

    **You can compare the Raw and Trimmed quality:**
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
   - Alternative tools: NanoPlot, pycoQC, PoreChop, FiltLong
   - Purpose: Evaluate read length distribution, quality scores, and sequencing statistics.

        ```bash
        fastplong -i <reads.fastq> -o <OUTPUT_FILENAME> -h qc.html -j qc.json -l 320
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

## Taxonomic Classification with Kraken2

### Goal:
Taxonomic classification of DNA or RNA sequences used to find exact matches in pre-built database of reference genomes.

Kraken2 using exact k-mer matches to compare sequences against a reference database, assigning them to specific taxonomic groups with high accuracy. Kraken2 is particularly useful for metagenomic analysis, enabling the identification and quantification of microbial communities in complex samples. It is widely used in clinical diagnostics for pathogen detection and in environmental studies for microbiome profiling. Its flexibility, speed, and ability to handle both short and long reads make it a valuable tool for researchers.


### **1. Installation**
- Follow the instructions to install Kraken2 from [GitHub](https://github.com/DerrickWood/kraken2)

- See [more details](https://ccb.jhu.edu/software/kraken2/) and publications about Kraken2.
- Kraken2 needs to be installed in its own repository folder, then scripts can be copied to another location if preferred.
- Scripts: `kraken2`, `kraken2-build`, `kraken2-inspect`

### **2. Build the Database**
This guide uses SILVA-138 SSURef-NR99 database, but Kraken2 has several other built-in reference and able to create from [external sources](https://benlangmead.github.io/aws-indexes/k2).

- **Create a database directory:**
   ```bash
   mkdir /kraken2_silva
   ```
- **Add SILVA sequences to the Kraken2 database:**
   ```bash
   kraken2-build --special silva --db /kraken2_silva
   ```

   - **Note:** This process may be **memory (RAM) intensive** and sufficient disk space necessary, due to the database size.

### **3. Run Kraken2 for Taxonomic Classification**


- **Paired-end reads (Short):**
  ```bash
  kraken2 --db </path/to/kraken2_silva> \
    --paired <trimmed_R1.fastq> <trimmed_R2.fastq> \
    --threads <CPU> \
    --confidence 0.2 \
    --output kraken2_output.txt \
    --report kraken2_report.tabular
  ```
- **Single-end reads (Long):**
  ```bash
  kraken2 --db </path/to/kraken2_silva> \
    <trimmed_reads.fastq> \
    --threads <CPU> \
    --confidence 0.2 \
    --output kraken2_output.txt \
    --report kraken2_report.tabular
  ```
    
    **Parameters:**
    
  `--db`: Built Kraken2 database location
    
  `--paired`: Paired-end reads, usually short-read sequences are paired-end
 
  `--threads`: Your CPU threads, if not sure, set to `8`
  
  `--confidence`: Confidence score threshold must be between 0 and 1, see more details below
    
  `--output`: Read sequence output

  `--report`: Taxonomic result, file name can be anything readable (.txt, etc.)

  
- **Confidence score**

    Confidence score is a numerical value that represents the reliability or certainty of a prediction, classification, or analysis result. It provides a measure of how likely a computational result is to be correct, helping researchers prioritize high-confidence findings and filter out uncertain ones.
  
    - Articles:
  
      - [Impact of database choice and confidence score on the performance of taxonomic classification using Kraken2](https://link.springer.com/article/10.1007/s42994-024-00178-0)
      - [An Alignment Confidence Score Capturing Robustness to Guide Tree Uncertainty](https://pmc.ncbi.nlm.nih.gov/articles/PMC2908709/)
      - [The standardisation of the approach to metagenomic human gut analysis: from sample collection to microbiome profiling](https://www.nature.com/articles/s41598-022-12037-3)

    
**IMPORTANT:** If you encounter a result like `100% classification`, it is a **false positive**, which is likely caused by incorrect database.

---

## **Refined Abundance Estimation with Bracken**

### Goal:

Bracken is used for refining taxonomic abundance estimates in metagenomic studies. It provides more accurate abundance estimates compared to raw Kraken2 outputs by resolving ambiguous classifications and redistributing reads to the most likely taxon. Bracken is computationally efficient and can handle large metagenomic datasets, making it suitable for high-throughput studies.

### **1. Get average read length**
- Tool: **[SeqKit](https://github.com/shenwei356/seqkit)**
- Alternative tools: NanoPlot, seqtk
- Purpose: Evaluate average read length for Bracken


- After installing SeqKit, navigate to the reads folder and run the following:
    ```bash
    seqkit stats *.fastq.gz > read_length_summary.txt
    ```
  - Your file extension might be different than `.fastq.gz`
- Once the script finished, open `read_length_summary.txt` and you can see the statistics  of each sample. As we need a single overall number for Bracken, need to process the `num_seqs` and `sum_len`.
- Formula:
    ```text
    Sum of "sum_len" / Sum of "num_seqs"
    ```
  - Use **Excel/Libre** if you need a fast solution, without writing a script.
  - This parameter **does not need to be strictly accurate** (feel free to round up to a whole number), Bracken needs only a ballpark-figure.

### **2. Install Bracken**
- [Download from GitHub](https://github.com/jenniferlu717/Bracken) and follow the documentation.

### **3. Run Bracken build**
- Purpose: To process the Kraken2 database to create files that describe the k-mer distribution for each taxonomic level (e.g., species, genus, family).
    
    ```bash
    bracken-build -l <READ_LEN> -d </path/to/kraken2_silva> -x </kraken2/script> -t <CPU>
    ```
    **Parameters:**
    
    `-l`: Average read length from SeqKit
    
    `-d`: Path to Kraken2 database

    `-x`: Kraken2 script file

    `-t`: Your CPU threads, if not sure, set to `8`

    - Optionally `-k` flag can be used to set the K-mer length that used to build the kraken database, if not set either, it uses `35` as default in Kraken2 and Bracken.
    - In practice, several builds can be created with different K-mer options for various project requirements, for example: `15, 30, 50, 75, 100, etc`.


### **4. Run Bracken**
- Generate the bracken report files for visualizations.
  
    ```bash
  bracken -d </path/to/kraken2_silva> \ 
  -i <kraken2_report_file> \ 
  -o bracken_report.bracken \ 
  -t <THRESHOLD> \ 
  -r <AVERAGE_READ_LEN> \ 
  -l <LEVEL>
   ```
  
    **Parameters:**

    `-d`: Path to Kraken2 database

    `-i`: Kraken2 report file (.tabular)

    `-o`: Output file name

    `-r`: Average read length from SeqKit, it **has to be** the same number as `bracken-build -l` parameter

    `-l`: Level of estimate abundance (Genus, Family etc.). Options: `D,P,C,O,F,G,S,S1`

    `-t`: Threshold to report only taxa that have an abundance of at least the set number percentage in the sample being analyzed
----


## **Visualization**

### Goal:

Visualizing the results from a Bracken report is essential for interpreting the taxonomic composition and abundance of species or taxa in a metagenomic sample.

Now the issue with visualization is, it is nearly impossible to give a general guideline, as it is always tailor-made for a specific project, based on the samples and requirements. However, here are some information that can help to get started.


### 1. Data Preparation
Before you attempt to visualize your results, make sure the generated Bracken file is in the **right format** for the selected tool. Some tools do that automatically, but many of them expect a certain input format.

- **Combine report files:**
  - This step is usually mandatory as you most likely have multiple samples, therefore multiple report files from each of them. Those files need to be merged into a single one in order to work with any tool.
  - Visualization tools often give information about the expected format, but generally a comma separated (.csv) text file is the common format.

- **Normalize data:**
   - It is a good practice to filter low-abundance taxa below a certain threshold, in order to work with a smaller data source, which is computationally less intensive.

Once the combined and optionally a filtered Bracken report file is ready, the visualization should be quick and ready to tweak for the requirements of the project. 

### 2. Tools

#### **Bash:**
While Bash itself is primarily a command-line tool and does not include direct visualization capabilities, you can use it to call other tools or manipulate files before visualizing them. A common approach is utilizing command-line utilities and file converters.

- **awk:** Can be used to filter and summarize the output data from Bracken reports.
- **sed:** Useful for text manipulation and data formatting in the report files.
- **gnuplot:** A powerful command-line plotting utility that can be scripted to visualize data, though it requires some setup to parse the specific format of Bracken output files.


#### **Command Line:**

- **Krona:** Interactive hierarchical pie charts for taxonomic data.

- **Pavian:** A web-based tool for interactive visualization and analysis of metagenomic classification results (supports Kraken/Bracken reports).

- **GraPhlAn:** A tool for creating circular taxonomic trees with abundance annotations.


#### **Python:**
Python is preferred approach to build entire pipelines and integrate the data manipulation and visualization part into.

- **Pandas:** For data manipulation of the Bracken report files, allowing for easy filtering and grouping.
- **Matplotlib:** A widely used plotting library for creating static, interactive, and animated visualizations in Python.
- **Seaborn:** A statistical data visualization library based on Matplotlib that provides a high-level interface for drawing attractive graphics.
- **Plotly:** An interactive graphing library that can create dashboards and dynamic visualizations.
- **Biopython:** If you need to handle biological sequence data along with Bracken visualization tasks.


#### **R:**
R is well-known for statistical analysis and comes with excellent visualization libraries:

- **ggplot2:** A popular R package for creating complex multi-layered graphics easily from data frames.
- **dplyr:** Used for data manipulation, often in conjunction with ggplot2 for plotting.
- **phyloseq:** Specifically designed for handling and visualizing microbiome data, which can incorporate Bracken results for taxonomic classification.


#### **Conda suite:**
Anaconda or BioConda is a distribution of Python and R programming languages designed for scientific computing, data science, machine learning and bioinformatics (BioConda). It simplifies package management and deployment. It is a more approachable solution than CLI tools, as all the necessary software are in the same place and provides a graphical interface.


#### **Online Tools:**
There are some online tools that can help visualize taxonomic relationships and abundance data, often with user-friendly interfaces.

- **Qiime2:** While primarily focused on microbial community analysis, Qiime2 has built-in visualization tools that can process outputs from tools like Bracken, provided that the data is formatted correctly.
- **MicrobiomeAnalyst:** A web-based platform for comprehensive microbiome data analysis and visualization. Simply upload the Bracken report file to the website and explore bar plots, heatmaps, alpha/beta diversity, etc.
- **Galaxy:** An open-source platform that allows users to run bioinformatics tools through a web interface. You can integrate Bracken with Galaxy for visualizations.
- **GraphPad Prism:** While primarily a statistical analysis software, it can be used interactively for data visualization if you input the Bracken output manually.


### **3. Common Plots**
- **Bar Plots:**
  - Show the relative abundance of taxa (e.g., species, genus, phylum) in a sample or across multiple samples.
  - Useful for comparing taxonomic composition between samples.

- **Heatmaps:**
  - Display the abundance of taxa across multiple samples, with colors representing abundance levels.
  - Useful for identifying patterns or clusters in large datasets.

- **Pie Charts:**
  - Represent the proportion of different taxa in a single sample.
  - Useful for a quick overview of dominant taxa.

- **Taxonomic Trees:**
  - Visualize the hierarchical relationships between taxa (e.g., phylum → genus → species) and their abundances.
  - Useful for exploring the structure of the microbial community.

- **Alpha Diversity Plots:**
  - Show the diversity within a sample (e.g., Shannon index, Simpson index).
  - Useful for comparing diversity across samples.

- **Beta Diversity Plots:**
  - Show differences in taxonomic composition between samples (e.g., PCA, PCoA, NMDS).
  - Useful for identifying sample groupings or outliers.

---

### Publications

[Metagenome analysis using the Kraken software suite](https://www.nature.com/articles/s41596-022-00738-y)

[Bracken: Estimating species abundance in metagenomics data](https://peerj.com/articles/cs-104/)

[Ultrafast and accurate 16S rRNA microbial community analysis using Kraken 2](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00900-2)

[Kraken: Ultrafast metagenomic sequence classification using exact alignments](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r46)

[Benchmarking taxonomic assignments based on 16S rRNA gene profiling](https://academic.oup.com/gigascience/article/7/5/giy054/4995265)

[Metagenomic profiling pipelines improve taxonomic classification for 16S amplicon sequencing data](https://www.nature.com/articles/s41598-023-40799-x)

[Studying long 16S rDNA sequences with ultrafast-metagenomic sequence classification using exact alignments (Kraken)](https://www.sciencedirect.com/science/article/abs/pii/S0167701216300112?via%3Dihub)

[Profiling of Oral Bacterial Communities](https://journals.sagepub.com/doi/10.1177/0022034520914594)

[Metagenomic evidence for a polymicrobial signature of sepsis](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000642)

[A microbial signature following bariatric surgery is robustly consistent across multiple cohorts](https://www.journalofinfection.com/article/S0163-4453(20)30287-5/fulltext)

