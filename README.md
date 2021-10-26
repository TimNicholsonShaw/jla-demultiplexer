# jla-demultiplexer: Internal JLA lab tool for processing gene-specific sequencing experiments

## Installation

```bash
pip install jla-demultiplexer
```

## Usage

fastqbreakdown: Removes duplicate reads and trims barcodes and randommer addition
```bash
fastqbreakdown -r1 [read 1 location] -r2 [read 2 location] \
  -r [length of random-mer] -b [barcode]
```

## Notes

- Random-mer length supplied should be 10 or 11
- barcode supplied should be barcode column + gene specific column from manifest

## Usage

demultiplexer: Use JLA standard manifests to align and create tail files. Very specific to our gene-specific sequencing methodology
```bash
demultiplexer -r1 [read 1 location] -r2 [read 2 location] \
  -m [manifest location] -e [ensids] -f [reference fasta]
```
## Notes

- Must provide either -e (comma separated) or -f, but not both
- outputs to location of read 1
- tail file output is compatible with Tailer-Analysis which can be found at https://timnicholsonshaw.shinyapps.io/tailer-analysis/

## Usage

jla-trim: Removes duplicate reads and trims randommer addition only. FASTQ output
```bash
jla-trim -r1 [read 1 location] -r2 [read 2 location] \
  -r [length of random-mer]
```

## Notes

- Random-mer length supplied should be 10 or 11



