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
bsub -q short -n 4 -W 2:00 -R "rusage[mem=8000]" \
>   scripts/align_trim_ampliconEZ.sh \
>   /home/cesar.bautistasotelo-umw/Genomes/hg38.fa \
>   /home/cesar.bautistasotelo-umw/projects/MM_analysis/1_ampliconEZ/1_mosaic_LEMD3_WDR37/30-1201770023/00_fastq \
>   /home/cesar.bautistasotelo-umw/projects/MM_analysis/1_ampliconEZ/1_mosaic_LEMD3_WDR37/30-1201770023/03_repo-test \
>   config/primers_example.tsv \
>   both
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

## 📊 Example QC: `samtools flagstat` (Untrimmed vs Trimmed)

After running the pipeline, you can quickly assess alignment quality using:

```
samtools flagstat SAMPLE_untrimmed.bam
samtools.flagstat SAMPLE_trimmed.bam
```

Untrimmed BAM example:
```
692,671 + 0 in total (QC-passed reads + QC-failed reads)
683,504 + 0 primary
665,659 + 0 mapped (96.10%)
581,132 + 0 properly paired (85.02%)
2,768 + 0 singletons (0.40%)
9,167 + 0 supplementary
```

Trimmed BAM example:
```
692,561 + 0 in total (QC-passed reads + QC-failed reads)
683,504 + 0 primary
669,447 + 0 mapped (96.66%)
647,098 + 0 properly paired (94.67%)
2,524 + 0 singletons (0.37%)
9,057 + 0 supplementary
```

Interpretation
	•	Read count stays constant → trimming preserves all read pairs.
	•	Mapped % increases after removing primer sequence.
	•	Properly paired % improves substantially (85% → 95%), which is expected when primer sequence is clipped off and insert sizes become more consistent.
	•	Slight reduction in singletons and supplementary alignments indicates cleaner alignments.

This improvement is typical for PCR amplicon sequencing, where primer sequences remain at the ends of reads and can interfere with alignment unless trimmed.

---

## 🤝 Author

Developed at **UMass Chan Medical School** for internal and collaborative use.

---
