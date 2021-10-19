# jla-demultiplexer: Internal JLA lab tool for processing gene-specific sequencing experiments

## Installation

```bash
pip install jla-demultiplexer
```

## Usage

fastqbreakdown: Removes duplicate reads and trims barcodes and randommer addition
```bash
fastqbreakdown -r1 [read 1 location] -r2 [read 2 location] -r [length of random-mer] -b [barcode]
```

## Notes
- Random-mer length supplied should be 10 or 11
- barcode supplied should be barcode column + gene specific column from manifest

## Coming soon
- Handling manifests
- Aligning
- Creating tailer-analysis compatible files

