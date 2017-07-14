# FAS
Analysis of Fatty Acid Syntases

## Installation

Dependencies:
 - Conda

```bash
  git clone git@github.com:thiesgehrmann/FAS.git
  cd FAS
  yes | conda create --name snakemake snakemake
  source activate snakemake
  snakemake --use-conda --cores 20 all
```
