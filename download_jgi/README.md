# Download JGI data

In order to make use of these scripts to download the published genomes from the JGI databases, make sure that you add yur JGI username and password in the `download.sh` script.
The python scripts rely only on standard libraries.

You can run the full pipeline with:

```bash
  ./download_jgi.sh

```

The pipeline produces a file called data.json, which contains a description of the files downloaded for each genome.
