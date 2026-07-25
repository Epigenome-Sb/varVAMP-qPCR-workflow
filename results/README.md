# Results

The workflow creates one result subdirectory per project:

```text
results/<project_name>/
```

## Generated outputs

Depending on the selected workflow options and whether VarVAMP identifies valid qPCR candidates, the project result directory may contain:

* the final multiple sequence alignment;
* the conservative consensus sequence;
* per-position conservation metrics;
* a conservation summary;
* VarVAMP primer and probe tables;
* BED coordinate files;
* PDF visualizations;
* PNG visualizations;
* the VarVAMP log.

Typical files include:

```text
<project_name>_alignment.fasta
<project_name>_consensus.fasta
<project_name>_conservation_by_position.csv
<project_name>_conservation_summary.txt

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

Not every VarVAMP output is guaranteed to be generated. Output availability depends on whether candidate primer/probe systems satisfy the selected design parameters and VarVAMP filters.

## Parameter-dependent results

Workflow results depend directly on the scientific parameters selected during execution.

For reproducibility, record the following parameters for every analysis.

### CD-HIT-EST

When sequence dereplication is enabled, record:

* identity threshold (`-c`);
* clustering mode (`-g`);
* strand comparison mode (`-r`);
* automatically selected word size (`-n`);
* whether CD-HIT-EST dereplication was enabled or skipped.

Different clustering parameters may change the number and composition of sequence representatives passed to MAFFT.

### MAFFT

Record the selected multiple sequence alignment strategy, for example:

```text
Auto
FFT-NS-1
FFT-NS-2
FFT-NS-i
NW-NS-2
NW-NS-i
L-INS-i
G-INS-i
E-INS-i
NW-NS-PartTree-1
```

The selected alignment strategy may affect alignment quality, gap placement, conservation estimates, consensus generation, and downstream VarVAMP designs.

### Conservation analysis

Record:

* minimum occupancy threshold;
* minimum major-base frequency threshold.

These parameters determine which alignment positions are classified as highly conserved and which positions contribute a nucleotide to the conservative consensus.

### VarVAMP

Record:

* consensus threshold (`-t`);
* maximum primer ambiguity (`-a`);
* maximum probe ambiguity (`-pa`);
* number of threads used.

These parameters may directly affect which primers, probes, and qPCR designs are generated.

## Reproducibility metadata

When preserving or publishing workflow results, also record:

* project name;
* input dataset name or version;
* source database;
* sequence retrieval date;
* accession numbers when available;
* inclusion and exclusion criteria;
* software versions;
* workflow version;
* execution date;
* all selected scientific parameters;
* any manual filtering or post-processing performed after the workflow.

For analyses intended for publication, it is recommended to retain a complete record of the parameter summary displayed by the workflow before execution.

## Interpretation

Generated outputs are computational predictions.

They do not, on their own, demonstrate:

* analytical specificity;
* analytical sensitivity;
* diagnostic sensitivity;
* diagnostic specificity;
* primer efficiency;
* probe performance;
* absence of off-target amplification;
* experimental validity;
* clinical validity.

Candidate primers and probes should therefore undergo additional in silico evaluation and appropriate experimental validation before diagnostic or research use.

## Git tracking

Generated project result directories are normally excluded from Git through:

```text
results/*
!results/README.md
```

This prevents large or dataset-specific result files from being committed accidentally while keeping this documentation file under version control.

Results intended for publication or public reproducibility may instead be archived separately, for example in a release, supplementary dataset, or scientific data repository.
