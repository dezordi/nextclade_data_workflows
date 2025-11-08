# ZikaV Genotypes Dataset Workflow

This repository contains the source code and datasets used to implement the Zika virus (zikaV) genotypes system in NextClade. The workflow is based on the organizational structure initially developed for the Chikungunya Virus by V-GEN-Lab team, available [here](https://github.com/V-GEN-Lab/nextclade-datasets-workflow/tree/main/chikV).

## Overview

This repository contains the workflow for generating datasets from the representative trees for Zika virus genotypes. These datasets are designed for integration into the NextClade tool, facilitating the implementation and analysis of Zika virus genotypes.

## Repository Contents

- **datasets/**: Contains the curated ZikaV sequence datasets.
- **scripts/**: Includes the scripts used for data processing and analysis.
- **config/**: Maps various ZikaV genotypes to the correct genomic annotation file.
- **resources/**: Files necessary for constructing and configuring the datasets.

## Usage

```bash
snakemake --cores 2
```

## Maintainer

- [Filipe Zimmer Dezordi](https://github.com/dezordi)
