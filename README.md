# 🐙 Kraken2 Custom Metagenomics Pipeline

A fully modular, automated pipeline designed to build **multiple Kraken2 custom databases**, run batch metagenomic classification, and generate **stacked taxonomic Excel reports**.  
Ideal for projects where you need independent databases (Plants, Insects, Bacteria, Viruses, etc.) without mixing references or wasting disk space.

This pipeline features:

### 🔥 Key Advantages

| Feature | Description |
|--------|-------------|
| **Zero-Duplication Taxonomy** | Only **one central taxonomy folder (~60GB)** is stored — all databases reuse it via **symlinks** |
| **Modular FASTA References** | Each project has its own reference folder (`P`LANTS`, `INSECTS`, ...) |
| **Argument-Based Execution** | Every script receives the **DB name as parameter** |
| **Batch Classification** | Automatically processes **all samples** in `raw_fastq/` |
| **Final Excel Matrix** | Consolidates reports into a **single stacked table** |

---

## 🧠 Architecture Overview

### 1. Modular FASTA Reference System

Each database reads only from its own reference folder:

```

data/fasta_ref/PLANTS/   → Genomes for PLANTS DB
data/fasta_ref/INSECTS/ → Genomes for INSECTS DB

````

No cross-contamination — perfect for multiple independent projects.

---

### 2. The "Zero-Duplication" Taxonomy Strategy

Standard Kraken2 DBs store local `taxonomy/` copies (~60GB per DB).  
Building 5 DBs = **300GB+ wasted space**.

#### Our solution:

- **📦 Centralized storage**: one `data/taxonomy/` directory shared by all DBs
- **🔗 Symlink-based linking**: each DB folder contains only a link to the master taxonomy
- **💾 <1KB storage per DB**
- **⬆ Auto-migration tool**: detects old taxonomy and offers to relocate & link automatically

---

### 3. Parametric Execution

Run scripts by specifying the database name:

```bash
./01_build_db.sh PLANTS
./01_build_db.sh INSECTS   # Reuses same taxonomy automatically
````

No code editing. No duplicated folders.

---

## 🚀 Features

* 🔧 **Custom DB Builder**
* 🧬 Modular **FASTA reference input**
* 🧠 **Central taxonomy with symlink reuse**
* 📊 **Stacked Excel output table**
* 🌀 **Batch paired-end processing**
* 🧽 Automatic filename cleaning for sample names
* 🧱 Scales easily to 10+ DBs without size explosion

---

## 📂 Repository Structure

```text
kraken_pipeline/
├── data/
│   ├── raw_fastq/             # Paired-end reads: *_1.fastq.gz & *_2.fastq.gz
│   ├── taxonomy/              # MASTER NCBI TAXONOMY (shared via symlinks)
│   ├── fasta_ref/             # Reference genomes grouped by DB name
│   │   ├── PLANTS/
│   │   └── INSECTS/
│   └── dbs/                   # Built Kraken2 databases
│       ├── PLANTS/            # Contains symlink → ../taxonomy/
│       └── INSECTS/
│
├── results/
│   ├── reports/               # Kraken2 report outputs
│   │   ├── PLANTS/
│   │   └── INSECTS/
│   └── final_tables/          # Excel summary files
│
├── scripts/
│   ├── 01_build_db.sh         # Builds DB from fasta_ref/<DB_NAME>
│   ├── 02_run_kraken.sh       # Classifies all samples
│   └── 03_generate_table.py   # Creates stacked Excel taxonomic matrix
│
├── requirements.txt
├── LICENSE
└── README.md
```

---

## 📖 Usage Guide

### **Phase 1 — Build a Database**

Create reference folder and add any number of `.fasta` genomes:

```bash
mkdir -p data/fasta_ref/PLANTS
cp genomes/*.fasta data/fasta_ref/PLANTS/
```

Build DB:

```bash
cd scripts/
./01_build_db.sh PLANTS
```

Script will:

1. Check for taxonomy
2. Download or migrate if needed
3. Create symlink into `dbs/PLANTS`
4. Build Kraken2 DB

---

### **Phase 2 — Classification**

Place reads inside `data/raw_fastq/`:

Naming required:

```
sampleA_1.fastq.gz
sampleA_2.fastq.gz
```

Run:

```bash
cd scripts/
./02_run_kraken.sh PLANTS
```
> **⚠️ Important Note:** If you run the script without any arguments (e.g., just `./02_run_kraken.sh`), it will look for a default database named **`CUSTOM_DB`**. To use your specific databases (like `FISH` or `PLANTS`), you **must** provide the name as shown.

Reports saved to:

```
results/reports/PLANTS/
```

---

### **Phase 3 — Reporting**

```bash
cd scripts/
python3 03_generate_table.py PLANTS
```

Output:

```
results/final_tables/PLANTS_Taxonomy_Stacked.xlsx
```

Example structure:

| Rank   | TaxID | Name       | Sample_A | Sample_B | Sample_C |
| ------ | ----- | ---------- | -------- | -------- | -------- |
| Phylum | 1234  | Arthropoda | 1200     | 500      | —        |
| Class  | 5678  | Insecta    | 800      | 300      | 50       |

---

## 🛠 Requirements

### System

```bash
conda install -c bioconda kraken2
# or
sudo apt install kraken2
```

### Python

```bash
pip install -r requirements.txt
# includes pandas + openpyxl
```

---

## ⚙ Configuration (Optional)

Inside scripts you can edit:

```bash
THREADS=12
DB_PATH="../data/dbs/PLANTS"
INPUT_DIR="../data/raw_fastq"
```

---

## 📝 License

MIT — Free for commercial & academic use.

---

### If you want I can also generate:

✔ badges
✔ usage GIFs
✔ example workflow diagram
✔ citation template for publications

Just ask. 🚀

```

---

Si deseas, puedo generar versión **más corta**, **con iconografía visual**, o un **README premium estilo profesional GitHub con banners**.
```
