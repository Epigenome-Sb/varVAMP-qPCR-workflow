# Results

The workflow creates one subdirectory per project:

```text
results/<project_name>/
```

This directory may contain:

- the final multiple sequence alignment;
- the conservative consensus;
- per-position conservation metrics;
- a conservation summary;
- VarVAMP primer and probe tables;
- BED coordinate files;
- PDF and PNG visualizations;
- the VarVAMP log.

Generated results are computational predictions. They do not demonstrate assay
specificity, diagnostic performance, or experimental validity on their own.

When committing results, also record:

- the input dataset version;
- software versions;
- command-line parameters;
- execution date;
- any manual filtering performed after the workflow.
