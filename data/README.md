# Input data

Place nucleotide FASTA files in this directory, or provide another path when
running the workflow.

Example:

```text
data/my_sequences.fasta
```

Raw FASTA files are ignored by the default `.gitignore` to avoid accidentally
publishing large, restricted, confidential, or insufficiently documented data.

For a reproducible public repository, provide an `accession_numbers.tsv` file
with fields such as:

```text
accession    genotype_or_lineage    database    retrieval_date
```

Also document:

- the source database;
- the exact retrieval date;
- search terms or query;
- inclusion criteria;
- exclusion criteria;
- sequence region or genome coverage;
- treatment of duplicates and incomplete sequences;
- any metadata used for stratified validation.

Only redistribute sequence data when their source terms and your institutional
or collaborator agreements permit it.
