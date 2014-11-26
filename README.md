quant3p
=======

A set of scripts for 3' RNA-seq quantification

## Install

To isntall run `python setup.py install`

## Example

We bundled a small example dataset.

```bash
cd example
```

There are two problems of 3' RNA-seq we are dealling with:
* Bad 3' annotaion.
* Relatively high number of multimappers near the 3' end.

### Fixing annotation

Basically we find peaks in RNA-seq and try to associate it with a close transcript.

First, we call peaks for positive and negative strands using `macs2-stranded` for
all sample combined to determine posible positions of unannotated exons.

```bash
macs2-stranded -n example ./bam/*.bam
```

It will produce a file `example_peaks.narrowPeak`.

Second, we fix annotation with RNA-seq peaks up to 5Kbp downstream of transcript's 3' end.
We also put newly added exons in a sepearate file `mm10.slice.extns.gtf`.

```bash
gtf-extend -g mm10.slice.gtf -p example_peaks.narrowPeak -o mm10.slice.extended.gtf --extns-out mm10.slice.extns.gtf
```

These way we added 4 likely exons that were not annotated.

### Fixing multimappers

Main idea here is that we expect RNA-seq reads to be aligned to the genes.
So, we discard alignments of multimappers to the unannotated regions and 
fix multimapping annotation (`NH` SAM field) for such reads. Some of them
can be uniquely mapped to one gene.

```bash
mkdir bam_fixed
for f in bam/*.bam
do
    echo "$f"
    fix-mm -g mm10.slice.gtf "$f" -o "${f/bam/bam_fixed}"
done
```

## Dependencies

`macs2-stranded` depends on `samtools` and `macs2` executables.  `macs2` version should be >= 2.0.10

`gtf-extend.py` and `fix-mm.py` needs `HTSeq` python package installed.

