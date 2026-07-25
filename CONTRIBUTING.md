# Contributing

Contributions are welcome through GitHub issues and pull requests.

## Reporting a problem

Please include:

* the operating system;
* Python version;
* CD-HIT-EST version;
* MAFFT version;
* VarVAMP version;
* the complete command used;
* the scientific parameters selected;
* the full error message;
* a minimal non-confidential FASTA example when possible.

For interactive runs, please also report the selected values for:

* CD-HIT-EST identity threshold (`-c`);
* CD-HIT-EST clustering mode (`-g`);
* CD-HIT-EST strand comparison (`-r`);
* MAFFT alignment strategy;
* VarVAMP consensus threshold (`-t`);
* primer ambiguity (`-a`);
* probe ambiguity (`-pa`).

Do not upload confidential, identifiable, restricted, sensitive, or unpublished
sequence data to a public GitHub issue.

## Pull requests

Before submitting a pull request:

1. create a focused branch;
2. keep changes limited to one clear purpose;
3. preserve command-line compatibility when possible;
4. update the README and other documentation for user-visible changes;
5. verify Python syntax;
6. test the workflow on a small nucleotide FASTA dataset;
7. test interactive parameter validation;
8. test CD-HIT-EST with dereplication enabled and disabled;
9. test at least one CD-HIT-EST clustering mode;
10. test at least two MAFFT alignment strategies when MAFFT-related code is modified;
11. test VarVAMP parameter selection when VarVAMP-related code is modified;
12. test `--skip-varvamp`;
13. test at least one fully non-interactive command-line execution;
14. verify that generated result files are written to the expected directories.

## Scientific parameter changes

Changes affecting scientific parameters should be documented clearly.

Contributors should explain:

* which parameter was added or modified;
* which external software option it corresponds to;
* the accepted values;
* whether the change affects interactive execution, command-line execution, or both;
* whether existing analyses may produce different results after the change.

Scientific parameters should not be silently changed without corresponding
documentation.

## Code quality

Contributions should:

* remain compatible with nucleotide FASTA input;
* provide clear error messages for invalid input;
* avoid hard-coded dataset-specific paths;
* avoid introducing unnecessary dependencies;
* preserve reproducibility whenever possible.

Before committing Python changes, syntax can be checked with:

```bash
python -m py_compile varvamp_qpcr_workflow.py
```

The available command-line options can be checked with:

```bash
python varvamp_qpcr_workflow.py --help
```

## Data and privacy

Example datasets used for testing should be public, synthetic, or otherwise
authorized for redistribution.

Do not commit:

* confidential sequence data;
* personally identifiable data;
* restricted datasets;
* credentials or API keys;
* institution-specific private paths;
* large generated intermediate files.

## License

By contributing, you agree that your contribution may be distributed under the
repository license:

**GNU General Public License v3.0 or later (GPL-3.0-or-later).**
