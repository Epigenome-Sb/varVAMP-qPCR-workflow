# Generic VarVAMP qPCR Workflow

A reproducible command-line workflow for nucleotide sequence dereplication,
multiple sequence alignment, conservation analysis, and qPCR primer/probe design.

The workflow combines:

- **CD-HIT-EST** for nucleotide sequence dereplication;
- **MAFFT** for multiple sequence alignment;
- a Python-based per-position conservation analysis;
- **VarVAMP** for qPCR primer and probe design;
- generation and collection of tabular and graphical results.

## Scope

This repository is designed for nucleotide FASTA datasets, including viral genomes
or comparable genomic regions. It is not intended for protein sequences.

The workflow performs computational design only. Candidate primers and probes must
still be evaluated for taxonomic specificity, genotype or lineage coverage,
secondary structures, and experimental performance.

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
│   └── README.md
├── notebook/
│   └── README.md
└── results/
    └── README.md
```

The `work/` directory is generated automatically and is not tracked by Git.

## Requirements

The workflow requires:

- Python 3.9 or later;
- CD-HIT-EST;
- MAFFT;
- VarVAMP;
- Biopython;
- pandas;
- matplotlib;
- Pillow;
- PyMuPDF.

A tested Conda environment is provided in `environment.yml`.

## Installation with Conda

From the repository root:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate varvamp-qpcr-workflow
```

Activation is optional. The script can also be run directly with `conda run`.

## Input data

The input must be a nucleotide FASTA file containing at least one sequence.

Accepted nucleotide symbols include standard DNA/RNA characters and IUPAC
ambiguity codes. Protein FASTA files are rejected.

Place the input file in `data/`, or provide any valid path when prompted.

Example:

```text
data/my_sequences.fasta
```

The repository does not need to contain the raw FASTA dataset. For public
scientific projects, consider sharing accession numbers, retrieval dates,
inclusion criteria, and a script or instructions for rebuilding the dataset.

## Interactive execution

Run:

```bash
python varvamp_qpcr_workflow.py
```

The program asks:

```text
Name or path of the input FASTA file:
```

Enter a filename or path such as:

```text
data/my_sequences.fasta
```

The program first checks the current directory and then the `data/` directory
when only a filename is supplied.

## Execution without activating Conda

```bash
conda run --no-capture-output -n varvamp-qpcr-workflow \
  python varvamp_qpcr_workflow.py
```

## Non-interactive execution

```bash
python varvamp_qpcr_workflow.py \
  --input data/my_sequences.fasta
```

Example with explicit project and output directories:

```bash
python varvamp_qpcr_workflow.py \
  --input data/my_sequences.fasta \
  --project-name example_project \
  --workdir work/example_project \
  --results results/example_project
```

## Main parameters

| Option | Default | Description |
|---|---:|---|
| `--identity` | `0.95` | CD-HIT-EST identity threshold |
| `--min-occupancy` | `0.95` | Minimum column occupancy for conservation |
| `--min-major-frequency` | `0.95` | Minimum major-base frequency |
| `--varvamp-threshold` | `0.95` | VarVAMP consensus threshold |
| `--primer-ambiguity` | `2` | Maximum ambiguous bases in primers |
| `--probe-ambiguity` | `2` | Maximum ambiguous bases in the probe |
| `--threads` | `8` | Threads used by VarVAMP |
| `--skip-varvamp` | disabled | Run dereplication, alignment, and conservation only |

Display all options with:

```bash
python varvamp_qpcr_workflow.py --help
```

## Output organization

For an input named `my_sequences.fasta`, the default project name is
`my_sequences`.

Intermediate files are written to:

```text
work/my_sequences/
```

Final outputs are written to:

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
```

VarVAMP output availability depends on whether candidate qPCR designs pass its
filters.

## Conservation definitions

For each alignment position, the workflow calculates:

- A, C, G, and T counts;
- gap and ambiguity counts;
- occupancy;
- major base and major-base frequency;
- Shannon entropy;
- strict conservation;
- threshold-based high conservation.

By default, a position contributes a nucleotide to the conservative consensus
when:

1. occupancy is at least 95%; and
2. the major-base frequency is at least 95%.

Positions that do not satisfy both criteria are represented by `N`.

Coordinates reported by the workflow are alignment coordinates and are not
automatically equivalent to positions in a reference genome.

## Reproducibility notes

Before reporting results:

1. record the sequence database and retrieval date;
2. retain the accession list and inclusion/exclusion criteria;
3. record software versions;
4. run the workflow from a clean output directory;
5. evaluate coverage using the original dataset, not only cluster representatives;
6. assess candidates by genotype, lineage, or other relevant biological grouping;
7. perform independent in silico specificity checks;
8. validate selected assays experimentally.

## Third-party software

This repository does not redistribute CD-HIT-EST, MAFFT, or VarVAMP. They are
installed separately through the environment definition.

VarVAMP is developed by its original authors and distributed under the GNU
General Public License. Cite VarVAMP and the other external tools in any
scientific work that uses them.

## Citation

Citation metadata are provided in `CITATION.cff`. GitHub will display a
**Cite this repository** option when this file is present on the default branch.

## License

The workflow code in this repository is licensed under the
**GNU General Public License v3.0 or later**. See `LICENSE`.

Input datasets and generated scientific results may have separate terms,
depending on their original sources and institutional requirements.
