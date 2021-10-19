# jla-demultiplexer: Internal JLA lab tool for processing gene-specific sequencing experiments

## Installation

```bash
pip install jla-demultiplexer
```

## Usage

fastqbreakdown: Removes duplicate reads and trims barcodes and randommer addition
```bash
fastqbreakdown -r1 [read 1 location] -r2 [read 2 location] -r [length of random-mer 10 or 11] -b [barcode, barcode+gene specific sequence from manifest]
```

## Coming soon
- Handling manifests
- Aligning
- Creating tailer-analysis compatible files

