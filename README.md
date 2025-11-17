# ampliconEZ-pipeline
A lightweight, cluster-friendly pipeline for aligning amplicon sequencing data with optional primer trimming.  
Designed for the UMass Chan **SCI (LSF) cluster**, but works anywhere with `bwa`, `samtools`, and `cutadapt`.

---

## 🧬 Overview

This pipeline processes paired-end amplicon sequencing reads by:

1. **Aligning raw (untrimmed) reads** to a reference genome  
2. Optionally **trimming PCR primers** using a gene → primer lookup table  
3. Optionally **aligning trimmed reads**  
4. Automatically **skipping completed samples** based on existing BAM + BAI files

It is designed to be:

- **Robust** (safe to re-run)
- **Modular** (supports per-gene primer definitions)
- **Cluster-compatible** (LSF)
- **Simple** (one Bash script + one TSV file)

---

## 📁 Directory Structure

```
ampliconEZ-pipeline/
├── README.md
├── LICENSE
├── .gitignore
├── config/
│   └── primers_example.tsv
└── scripts/
    └── align_trim_ampliconEZ.sh
```

---

## 🧪 Primer Configuration

Primer trimming uses a **TSV file** with the following format:

```
gene    fwd_primer                  rev_primer
WDR37   AAAACGGCACATGATCCCCT        TTCCAAATCTCTCACCACTTCT
LEMD2   GGTCTAGCGTCTAGCTAGCTA       TCCGATCGATCGACTGACTGA
```

### Notes:
- The **gene name must appear in the FASTQ filename**  
  Example:  
  `MMgDNA-WDR37_S1_L001_R1_001.fastq.gz` → matches gene `WDR37`
- The script automatically selects the correct primer pair.

---

## 🎛 Modes of Operation

| TRIM_MODE | Behavior |
|----------|----------|
| `none` *(default)* | Only align untrimmed reads |
| `trimmed` | Trim primers + align trimmed reads only |
| `both` | Align untrimmed AND trimmed reads |

Primer trimming modes require a TSV.

---

## 🔧 Requirements

Modules needed:

```
module load bwa/0.7.17
module load samtools/1.3
module load cutadapt/4.1
```

Reference must be indexed:

```
bwa index hg38.fa
samtools faidx hg38.fa
```

---

## 🚀 Running on SCI Cluster (LSF)

Submit the script:

```bash
bsub < scripts/align_trim_ampliconEZ.sh
```

You can either:

- Edit defaults inside the script, OR  
- Pass arguments manually (see below)

---

## 🖥 Manual Run (no LSF)

```bash
bash scripts/align_trim_ampliconEZ.sh     /path/to/hg38.fa     /path/to/fastqs     /path/to/output     config/primers_example.tsv     both
```

### Arguments:

| Position | Description |
|----------|-------------|
| 1 | Reference FASTA |
| 2 | FASTQ directory |
| 3 | Output directory |
| 4 | Primer TSV (required for trimmed modes) |
| 5 | `none` \| `trimmed` \| `both` |

---

## 📤 Output Files

If untrimmed:

```
SAMPLE_untrimmed.bam
SAMPLE_untrimmed.bam.bai
```

If trimmed:

```
trimmed/SAMPLE_R1_trimmed.fastq.gz
trimmed/SAMPLE_R2_trimmed.fastq.gz
trimmed/SAMPLE_cutadapt.log
SAMPLE_trimmed.bam
SAMPLE_trimmed.bam.bai
```

---

## 🔁 Resume-Safe Execution

The pipeline skips any existing BAM + BAI pairs:

```
[SKIP] SAMPLE trimmed alignment already exists.
```

Safe to re-run after failures or interruptions.

---

## 🤝 Author

Developed at **UMass Chan Medical School** for internal and collaborative use.

---
