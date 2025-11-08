## Create Zika test dataset

Scripts and files to create a test dataset for Zika virus

### References

The `scripts/generate_from_genbank.py` was modified from [nextstrain/nextclade_data/docs/example-workflow/scripts/generate_from_genbank.py](https://github.com/nextstrain/nextclade_data/blob/master/docs/example-workflow/scripts/generate_from_genbank.py).


The `get_test_dataset.smk` was modified from [nextstrain/nextclade_data/docs/example-workflow/Snakefile](https://github.com/nextstrain/nextclade_data/blob/master/docs/example-workflow/Snakefile).

### Usage

```bash
snakemake -s get_test_dataset.smk --cores 2
```

When prompted interactively, enter `2` and then `product` to instruct the `generate_from_genbank.py` script to split the polyprotein CDS into `mature_protein` regions.

### Outputs

The workflow will produce a directory structure like this:

```bash
├── data                              # Sequences and metadata from Zika virus genomes present on NCBI
│   ├── metadata.tsv
│   └── sequences.fasta
├── dataset_test                      # Sequence, genbank file and gff file from Reference genome 
│   ├── annotation.gff
│   ├── genome_annotation.gff3
│   ├── reference.fasta
│   └── reference.gb
├── results_dataset_test              # Files produced during the dataset creation (mainly augur outputs)
│   ├── aligned.fasta
│   ├── auspice.json
│   ├── branch_lengths.json
│   ├── example_sequences.fasta
│   ├── ...
├── test_out                          # Nextclade outputs from a test with nextclade run with the new dataset
│   ├── nextclade.aligned.fasta
│   ├── nextclade.auspice.json
│   ├── ...
├── zika_dataset_test                  # The dataset itself
│   ├── CHANGELOG.md
│   ├── genome_annotation.gff3
│   ├── pathogen.json
│   ├── README.md
│   ├── reference.fasta
│   ├── sequences.fasta
│   └── tree.json
```

Some of the files produced here can be used to configure the official dataset:

- data/metadata.tsv                             -> Remove from resources_test/exclude.txt and cut first 4 columns to produce `nextclade_data_workflows/zikaV/resources/raw/metadata.tsv`
- data/sequences.fasta                          -> Remove from resources_test/exclude.txt to produce `nextclade_data_workflows/zikaV/resources/raw/sequence.fasta`
- dataset_test/reference.fasta                  -> `nextclade_data_workflows/zikaV/resources/zikaV/reference.fasta`
- dataset_test/genome_annotation.gff3           -> `nextclade_data_workflows/zikaV/resources/zikaV/annotation.gff3`
- dataset_test/reference.gb                     -> `nextclade_data_workflows/zikaV/resources/reference.gb`
- zika_dataset_test/tree.json                   -> Open this file in Auspice (e.g., run `auspice view --datasetDir ./zika_dataset_test`) and inspect the tree to identify unique mutations for each clade; use these mutations to populate `nextclade_data_workflows/zikaV/resources/clades.tsv`. For more details, see the [Auspice documentation](https://docs.nextstrain.org/projects/auspice/en/latest/).
- zika_dataset_test/tree.json                   -> Check with auspice and get unique mutations to configure `nextclade_data_workflows/zikaV/resources/clades.tsv`
