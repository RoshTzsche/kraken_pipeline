# 🐙 Kraken2 Custom Metagenomics Pipeline

A modular, automated pipeline designed to build custom **Kraken2** databases, perform batch classification of Paired-End (PE) FASTQ sequences, and generate consolidated taxonomic reports in a clean Excel format.

This repository is structured to handle large datasets efficiently, keeping raw data, scripts, and results organized while minimizing disk usage.

---

## 🚀 Features

- **Custom DB Builder:** Compiles local `.fasta` files into a Kraken2 library. Includes an **interactive check** to automatically download NCBI taxonomy if it's missing.
- **Batch Classification:** Automatically processes all paired-end samples in the input directory.
- **Storage Efficiency:** Generates lightweight Kraken2 reports (`_report.txt`) without storing massive, unnecessary read-by-read output files.
- **Stacked Taxonomy Output:** Parses individual sample reports into a single, consolidated **Stacked Excel Matrix** containing Ranks, TaxIDs, and read counts for every sample.

---

## 📂 Repository Structure

The project uses a strict directory structure to ensure scripts run correctly.  
**Note:** Large data files (FASTQ, DBs) are excluded from the repository via `.gitignore`.

```text
kraken_pipeline/
├── data/
│   ├── raw_fastq/          # PLACE YOUR INPUT FILES HERE (*_1.fastq.gz, *_2.fastq.gz)
│   ├── fasta_ref/          # Place reference .fasta genomes here for building DB
│   └── dbs/                # Destination for the built Kraken2 database
│
├── results/
│   ├── reports/            # Generated text reports from Kraken2
│   └── final_tables/       # Final Excel summary files
│
├── scripts/
│   ├── 01_build_db.sh       # Builds the database (Interactive)
│   ├── 02_run_kraken.sh     # Runs classification on all samples
│   └── 03_generate_table.py # Compiles reports into a single Excel table
│
├── requirements.txt         # Python dependencies
├── .gitignore               # Ignores heavy data files
├── LICENSE                  # MIT License
└── README.md                # Documentation
````

---

## 🛠️ Prerequisites

### 1. System Requirements

* Linux/Unix environment (or WSL on Windows)
* Kraken2 installed and added to system PATH
  Install using:

```bash
conda install -c bioconda kraken2
# o
sudo apt install kraken2
```

### 2. Python Dependencies

The reporting script requires Python 3 with:

```bash
pip install -r requirements.txt
# Or manually:
pip install pandas openpyxl
```

---

## 📖 Usage Instructions

### **Step 1: Build the Database**

> Skip if you already have a compiled database.

1. Place `.fasta / .fa` reference genomes in `data/fasta_ref/`
2. Run:

```bash
cd scripts/
./01_build_db.sh
```

The script will check for NCBI taxonomy:

| Option | Action                                                        |
| ------ | ------------------------------------------------------------- |
| **y**  | Downloads automatically (requires internet)                   |
| **n**  | Stop & manually provide taxonomy inside `data/dbs/CUSTOM_DB/` |

---

### **Step 2: Run Classification**

1. Place paired-end reads in `data/raw_fastq/`

   * Naming required: `*_1.fastq.gz` & `*_2.fastq.gz`

2. Run:

```bash
cd scripts/
./02_run_kraken.sh
```

Reports will be saved in `results/reports/`
Sample names are parsed automatically (cleaning prefixes/suffixes).

---

### **Step 3: Generate Summary Table**

```bash
cd scripts/
python3 03_generate_table.py
```

### 📊 Output Format

Final table: `results/final_tables/Taxonomy_Stacked_Results.xlsx`

| Rank   | TaxID | Name       | Sample_A | Sample_B | Sample_C |
| ------ | ----- | ---------- | -------- | -------- | -------- |
| Phylum | 1234  | Arthropoda | 500      | 1200     | —        |
| Class  | 5678  | Insecta    | 450      | 100      | 50       |
| Order  | 9101  | Coleoptera | 300      | 80       | 10       |

---

## ⚙️ Configuration

Modify paths and settings inside scripts, example inside `02_run_kraken.sh`:

```bash
DB_PATH="../data/dbs/CUSTOM_DB"
THREADS=12
INPUT_DIR="../data/raw_fastq"
```

---

## 📝 License

This project is licensed under the **MIT License** — see `LICENSE` for details.

---