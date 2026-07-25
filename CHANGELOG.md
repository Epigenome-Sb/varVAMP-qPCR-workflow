# Changelog

All notable changes to this project will be documented in this file.

## [1.1.0] - 2026-07-25

### Added

* Interactive CD-HIT-EST parameter checkpoints.
* Optional CD-HIT-EST dereplication.
* Explicit CD-HIT-EST identity threshold selection.
* CD-HIT-EST fast (`-g 0`) and accurate (`-g 1`) clustering modes.
* CD-HIT-EST strand comparison selection using `-r`.
* Automatic CD-HIT-EST word-size selection according to the chosen identity threshold.
* Interactive MAFFT alignment strategy selection.
* Support for the following MAFFT strategies:

  * Auto;
  * FFT-NS-1;
  * FFT-NS-2;
  * FFT-NS-i with 2 refinement cycles;
  * FFT-NS-i with up to 1000 refinement cycles;
  * NW-NS-2;
  * NW-NS-i with 2 refinement cycles;
  * NW-NS-i with up to 1000 refinement cycles;
  * L-INS-i;
  * G-INS-i;
  * E-INS-i;
  * NW-NS-PartTree-1.
* Warnings for computationally intensive MAFFT strategies when the dataset size is outside their typical use range.
* Interactive VarVAMP parameter checkpoints.
* Explicit VarVAMP consensus threshold (`-t`) selection.
* Explicit primer ambiguity (`-a`) selection.
* Explicit probe ambiguity (`-pa`) selection.
* Final scientific parameter summary before workflow execution.
* Explicit user confirmation before starting the analysis.
* Support for fully non-interactive execution through command-line arguments.

### Changed

* Removed fixed interactive defaults for the main CD-HIT-EST and VarVAMP scientific parameters.
* MAFFT is no longer restricted to the `--auto` strategy.
* CD-HIT-EST execution now exposes clustering mode and strand comparison options.
* The workflow now requires explicit scientific parameter selection during interactive execution.
* Documentation updated to reflect the new parameter checkpoints and MAFFT strategies.
* Repository structure updated to remove the notebook component.
* Conda environment simplified by removing JupyterLab and IPython kernel dependencies.

### Improved

* Reproducibility by requiring users to explicitly review the scientific parameters used for each analysis.
* Input validation for interactive numerical parameters.
* Error handling for invalid parameter selections.
* Transparency of CD-HIT-EST, MAFFT, and VarVAMP execution settings.

## [1.0.0] - 2026-07-22

### Added

* Generic nucleotide FASTA input selection.
* CD-HIT-EST sequence dereplication.
* MAFFT multiple sequence alignment.
* Per-position conservation analysis.
* Conservative consensus generation.
* VarVAMP qPCR primer and probe design.
* Result collection and graphical exports.
* English command-line interface and documentation.
