# Generic VarVAMP qPCR Workflow

A reproducible and interactive command-line workflow for nucleotide sequence dereplication, multiple sequence alignment, conservation analysis, and qPCR primer/probe design.

The workflow combines:

* **CD-HIT-EST** for nucleotide sequence dereplication;
* **MAFFT** for multiple sequence alignment;
* a Python-based per-position conservation analysis;
* **VarVAMP** for qPCR primer and probe design;
* interactive scientific parameter checkpoints;
* generation and collection of tabular and graphical results.

## Scope

This repository is designed for nucleotide FASTA datasets, including viral genomes, genes, genomic regions, or other comparable nucleotide sequence collections.

It is not intended for protein sequences.

The workflow performs computational assay design only. Candidate primers and probes should subsequently be evaluated for:

* taxonomic specificity;
* genotype, lineage, or strain coverage;
* primer/probe mismatches;
* secondary structures;
* potential off-target amplification;
* experimental qPCR performance.

## Workflow overview

```text
Input nucleotide FASTA
        │
        ▼
FASTA validation
        │
        ▼
CD-HIT-EST checkpoint
        │
        ├── Run dereplication? yes/no
        ├── Identity threshold (-c)
        ├── Clustering mode (-g)
        └── Strand comparison (-r)
        │
        ▼
MAFFT checkpoint
        │
        └── Alignment strategy selection
        │
        ▼
Multiple sequence alignment
        │
        ▼
Conservation analysis
        │
        ▼
VarVAMP checkpoint
        │
        ├── Consensus threshold (-t)
        ├── Primer ambiguity (-a)
        └── Probe ambiguity (-pa)
        │
        ▼
Scientific parameter summary
        │
        ▼
User confirmation
        │
        ▼
qPCR design and result generation
```

The workflow does not silently impose fixed scientific choices for CD-HIT-EST, MAFFT, or VarVAMP during interactive execution. The user explicitly selects the relevant parameters before analysis begins.

## Repository structure

```text
.
├── varvamp_qpcr_workflow.py
├── README.md
├── LICENSE
├── CITATION.cff
├── environment.yml
├── .gitignore
├── CONTRIBUTING.md
├── CHANGELOG.md
├── data/
│   ├── README.md
│   └── accession_numbers.tsv
└── results/
    └── README.md
```

The `work/` directory is generated automatically during execution and is not intended to be tracked by Git.

## Requirements

The workflow requires:

* Python 3.9 or later;
* CD-HIT-EST;
* MAFFT;
* VarVAMP;
* Biopython;
* pandas;
* matplotlib;
* Pillow;
* PyMuPDF.

A Conda environment definition is provided in `environment.yml`.

## Installation with Conda

From the repository root:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate varvamp-qpcr-workflow
```

Environment activation is optional. The workflow can also be executed directly with `conda run`.

## Input data

The input must be a nucleotide FASTA file containing at least one sequence.

Accepted characters include standard DNA/RNA nucleotides and IUPAC ambiguity codes.

Protein FASTA files are rejected.

The FASTA file may be stored in the `data/` directory:

```text
data/my_sequences.fasta
```

or supplied from any valid location on the computer.

The repository does not need to contain the raw FASTA dataset.

For reproducible public studies, users are encouraged to document:

* accession numbers;
* source database;
* retrieval date;
* search strategy;
* inclusion criteria;
* exclusion criteria;
* genotype, lineage, or strain information;
* sequence completeness;
* genomic region analysed.

## Interactive execution

Run:

```bash
python varvamp_qpcr_workflow.py
```

The workflow first asks for the input file:

```text
Name or path of the input FASTA file:
```

Example:

```text
data/my_sequences.fasta
```

When only a filename is provided, the workflow searches the current directory and then the `data/` directory.

After validating the FASTA file, the workflow guides the user through the scientific parameter checkpoints.

---

# CD-HIT-EST checkpoint

The workflow first asks whether sequence dereplication should be performed:

```text
Run CD-HIT-EST sequence dereplication? [y/n]:
```

If CD-HIT-EST is selected, the user must choose the following parameters.

## Identity threshold

```text
CD-HIT-EST identity threshold (-c) [0.75-1.0]:
```

The value determines the minimum sequence identity required for clustering.

Example:

```text
0.95
```

The workflow automatically selects a compatible CD-HIT-EST word size (`-n`) according to the selected identity threshold.

## Clustering mode

The workflow then asks:

```text
CD-HIT-EST clustering mode checkpoint
=====================================

1. Fast mode (-g 0)
2. Accurate mode (-g 1)
```

### Fast mode

```text
-g 0
```

A sequence is assigned to the first cluster representative satisfying the selected identity threshold.

This approach is computationally faster.

### Accurate mode

```text
-g 1
```

A sequence is compared against the available representatives and assigned to the most similar qualifying cluster.

This approach may provide more precise cluster assignment but requires additional computation.

## Strand comparison

The workflow also controls nucleotide strand comparison:

```text
CD-HIT-EST strand checkpoint
============================

1. Compare both +/+ and +/- orientations (-r 1)
2. Compare only the same orientation, +/+ (-r 0)
```

For datasets in which sequence orientation may vary, comparison of both strands may be appropriate.

If all sequences are already consistently oriented, same-strand comparison may be sufficient.

## Skipping dereplication

CD-HIT-EST can be skipped entirely.

In that case, all input sequences are passed directly to MAFFT.

---

# MAFFT alignment checkpoint

The workflow requires an explicit multiple sequence alignment strategy.

The available strategies are:

```text
1. Auto
2. FFT-NS-1
3. FFT-NS-2
4. FFT-NS-i (2 cycles)
5. FFT-NS-i (up to 1000 cycles)
6. NW-NS-2
7. NW-NS-i (2 cycles)
8. NW-NS-i (up to 1000 cycles)
9. L-INS-i
10. G-INS-i
11. E-INS-i
12. NW-NS-PartTree-1
```

## Auto

MAFFT automatically selects an appropriate strategy according to the dataset characteristics and size.

## FFT-NS-1

Very fast progressive alignment intended primarily for large datasets.

## FFT-NS-2

Fast progressive alignment using two guide-tree calculations.

## FFT-NS-i

FFT-based alignment followed by iterative refinement.

Two configurations are available:

```text
FFT-NS-i (2 cycles)
FFT-NS-i (up to 1000 cycles)
```

The 1000-cycle configuration provides more intensive refinement while remaining more scalable than the computationally expensive pairwise strategies.

## NW-NS-2

Progressive alignment without FFT approximation.

## NW-NS-i

Iterative refinement without FFT approximation.

Two configurations are available:

```text
NW-NS-i (2 cycles)
NW-NS-i (up to 1000 cycles)
```

## L-INS-i

```text
--localpair --maxiterate 1000
```

High-accuracy strategy based on local pairwise alignment followed by iterative refinement.

It is suitable when sequences share a locally alignable region but may contain variable flanking regions.

This strategy is computationally intensive and is generally intended for relatively small datasets.

## G-INS-i

```text
--globalpair --maxiterate 1000
```

High-accuracy strategy based on global pairwise alignment followed by iterative refinement.

It is appropriate when sequences:

* are homologous across most of their length;
* represent the same genomic region;
* have similar lengths;
* are expected to be globally alignable.

Because global pairwise calculations are computationally demanding, G-INS-i is primarily intended for relatively small datasets.

## E-INS-i

```text
--genafpair --maxiterate 1000
```

Designed for sequences containing several conserved regions separated by long variable or difficult-to-align regions.

## NW-NS-PartTree-1

A scalable strategy intended for extremely large datasets.

The workflow displays a warning when a selected MAFFT strategy is not normally appropriate for the size of the current dataset and requires confirmation before continuing.

---

# Conservation analysis

After multiple sequence alignment, the workflow performs a per-position conservation analysis.

For each alignment column, it calculates:

* A, C, G, and T counts;
* number of gaps;
* number of ambiguous characters;
* occupancy;
* major base;
* major-base frequency;
* Shannon entropy;
* strict conservation;
* threshold-based conservation.

## Occupancy threshold

The default conservation-analysis occupancy threshold is:

```text
0.95
```

A position must therefore be represented in at least 95% of the analysed sequences to satisfy this criterion.

This threshold can be changed using:

```bash
--min-occupancy
```

## Major-base frequency

The default major-base frequency threshold is:

```text
0.95
```

This can be changed using:

```bash
--min-major-frequency
```

A position contributes its major nucleotide to the conservative consensus when both criteria are satisfied.

Positions failing either criterion are represented by:

```text
N
```

The conservation-analysis thresholds are independent of the VarVAMP consensus threshold.

Coordinates generated during this analysis correspond to multiple-sequence-alignment coordinates and are not automatically equivalent to coordinates in a reference genome.

---

# VarVAMP checkpoint

Before qPCR design, the workflow explicitly requests the main VarVAMP parameters.

No interactive default values are imposed for these parameters.

## Consensus threshold

```text
VarVAMP consensus threshold (-t):
```

Example:

```text
0.95
```

Higher values impose a more stringent consensus requirement.

## Primer ambiguity

```text
Maximum ambiguous bases in each primer (-a):
```

Example:

```text
2
```

This controls the maximum number of IUPAC ambiguous positions permitted in each primer.

A value of:

```text
0
```

prohibits ambiguous bases in primers.

## Probe ambiguity

```text
Maximum ambiguous bases in the probe (-pa):
```

Example:

```text
2
```

This controls the maximum number of ambiguous positions permitted in the qPCR probe.

---

# Final parameter confirmation

Before running the main analysis, the workflow displays all selected scientific parameters.

Example:

```text
Scientific parameter summary
============================

CD-HIT-EST identity (-c): 0.95
CD-HIT-EST mode: accurate
CD-HIT-EST strand comparison: both
CD-HIT-EST word size (-n): 10

MAFFT strategy: fft-ns-i-1000

VarVAMP consensus threshold (-t): 0.95
Primer ambiguity (-a): 2
Probe ambiguity (-pa): 2
```

The user must explicitly confirm:

```text
Run the workflow with these parameters? [y/n]:
```

The analysis begins only after confirmation.

---

# Execution without activating Conda

The workflow can be executed without manually activating the environment:

```bash
conda run --no-capture-output -n varvamp-qpcr-workflow \
  python varvamp_qpcr_workflow.py
```

Interactive questions remain available when `--no-capture-output` is used.

---

# Non-interactive execution

All required scientific parameters can also be supplied directly on the command line.

Example:

```bash
python varvamp_qpcr_workflow.py \
  --input data/my_sequences.fasta \
  --identity 0.95 \
  --cdhit-mode accurate \
  --cdhit-strand both \
  --mafft-strategy fft-ns-i-1000 \
  --varvamp-threshold 0.95 \
  --primer-ambiguity 2 \
  --probe-ambiguity 2 \
  --threads 8
```

With explicit project directories:

```bash
python varvamp_qpcr_workflow.py \
  --input data/my_sequences.fasta \
  --project-name example_project \
  --workdir work/example_project \
  --results results/example_project \
  --identity 0.95 \
  --cdhit-mode accurate \
  --cdhit-strand both \
  --mafft-strategy fft-ns-i-1000 \
  --varvamp-threshold 0.95 \
  --primer-ambiguity 2 \
  --probe-ambiguity 2
```

To skip CD-HIT-EST:

```bash
python varvamp_qpcr_workflow.py \
  --input data/my_sequences.fasta \
  --skip-cdhit \
  --mafft-strategy auto \
  --varvamp-threshold 0.95 \
  --primer-ambiguity 2 \
  --probe-ambiguity 2
```

To skip VarVAMP and perform only preprocessing, alignment, and conservation analysis:

```bash
python varvamp_qpcr_workflow.py \
  --input data/my_sequences.fasta \
  --identity 0.95 \
  --cdhit-mode accurate \
  --cdhit-strand both \
  --mafft-strategy auto \
  --skip-varvamp
```

Display all available options with:

```bash
python varvamp_qpcr_workflow.py --help
```

---

# Main command-line parameters

| Option                  | Interactive default | Description                                 |
| ----------------------- | ------------------: | ------------------------------------------- |
| `--input`               |                none | Input nucleotide FASTA file                 |
| `--project-name`        |      FASTA filename | Project identifier                          |
| `--identity`            |                none | CD-HIT-EST identity threshold               |
| `--skip-cdhit`          |            disabled | Skip sequence dereplication                 |
| `--cdhit-mode`          |                none | `fast` or `accurate` clustering             |
| `--cdhit-strand`        |                none | `both` or `same` strand comparison          |
| `--mafft-strategy`      |                none | MAFFT alignment strategy                    |
| `--min-occupancy`       |              `0.95` | Conservation-analysis occupancy threshold   |
| `--min-major-frequency` |              `0.95` | Conservation major-base frequency threshold |
| `--varvamp-threshold`   |                none | VarVAMP consensus threshold                 |
| `--primer-ambiguity`    |                none | Maximum ambiguous bases per primer          |
| `--probe-ambiguity`     |                none | Maximum ambiguous bases in the probe        |
| `--threads`             |                 `8` | Number of computational threads             |
| `--skip-varvamp`        |            disabled | Skip VarVAMP qPCR design                    |

---

# Output organization

For an input named:

```text
my_sequences.fasta
```

the automatically generated project name is:

```text
my_sequences
```

Intermediate files are written to:

```text
work/my_sequences/
```

Final results are written to:

```text
results/my_sequences/
```

Typical outputs include:

```text
my_sequences_alignment.fasta
my_sequences_consensus.fasta
my_sequences_conservation_by_position.csv
my_sequences_conservation_summary.txt

qpcr_primers.tsv
qpcr_design.tsv
primers.bed
amplicons.bed

amplicon_plot.pdf
per_base_mismatches.pdf

amplicons_overview.png
primers_overview.png

ambiguous_consensus.fasta
varvamp_log.txt
```

VarVAMP output availability depends on whether candidate qPCR systems satisfy the selected filters.

---

# Reproducibility

Because CD-HIT-EST, MAFFT, and VarVAMP parameters can now be selected interactively, the exact parameter configuration used for an analysis should always be recorded.

Before reporting or publishing results:

1. record the sequence database and retrieval date;
2. retain accession numbers and inclusion/exclusion criteria;
3. record software versions;
4. record the CD-HIT-EST identity threshold;
5. record the CD-HIT-EST clustering mode;
6. record the CD-HIT-EST strand-comparison mode;
7. record the MAFFT strategy;
8. record the conservation-analysis thresholds;
9. record the VarVAMP consensus threshold;
10. record primer and probe ambiguity limits;
11. record the number of threads;
12. run the workflow from a clean output directory;
13. evaluate assay coverage using the original dataset, not only cluster representatives;
14. evaluate candidates by genotype, lineage, strain, or other biologically relevant grouping;
15. perform independent in silico specificity analysis;
16. experimentally validate selected assays.

Computational predictions alone do not establish diagnostic sensitivity, specificity, amplification efficiency, or clinical validity.

---

# Third-party software

This repository does not redistribute CD-HIT-EST, MAFFT, or VarVAMP.

These programs are installed separately through the provided environment definition.

Users should cite the original software publications when reporting analyses performed with:

* CD-HIT/CD-HIT-EST;
* MAFFT;
* VarVAMP.

VarVAMP is developed and maintained by its original authors and is distributed under the GNU General Public License.

---

# Citation

Citation metadata for this workflow are provided in:

```text
CITATION.cff
```

When this file is present on the default GitHub branch, GitHub can display a:

```text
Cite this repository
```

option.

---

# License

The workflow code in this repository is licensed under the:

**GNU General Public License v3.0 or later**

See:

```text
LICENSE
```

Input sequence datasets and generated scientific results may be subject to separate terms depending on their source databases, institutional policies, collaborators, or data-use agreements.
