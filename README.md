quant3p
=======

A set of scripts for 3' RNA-seq quantification

## Install

To install run `python setup.py install`.

To install locally, without root privelegies, run `python setup.py install --user`. Don't forget to add `.local/bin` to your `PATH`. You can do it by adding line `export PATH="$HOME/.local/bin:$PATH"` to the `~/.profile` file.

## Quick start

You can quickly run analysis for the example data.

```bash
quant3p -n example -g example/mm10.slice.gtf example/bam/*.bam 
```

Parameters are:
* `-n NAME` sets the name of experiment to be `NAME`,
* `-g GTF` tells to use `GTF` as the annotation,
* `example/bam/*.bam` is a list of bam-files to process.

This command produces one file `example.cnt` that contatins gene counts, that looks like this:
```
                        2h_rep1  2h_rep2  4h_rep1  4h_rep2
100039246               0        0        0        0
11744                   112      118      72       86
14673                   105      114      52       67
211623                  0        0        0        0
231842                  47       39       16       17
232157                  1799     1947     2060     2075
66039                   1        2        0        0
78653                   147      148      90       131
NM_001270503_dup1       0        0        0        0
NM_001270503_dup2       0        0        0        0
NM_207229_dup1          0        0        0        0
NM_207229_dup2          0        0        0        0
__no_feature            463      489      389      433
__ambiguous             0        0        0        0
__too_low_aQual         0        0        0        0
__not_aligned           0        0        0        0
__alignment_not_unique  289      190      231      194
```

If you want to look at intermediate files, you can use parameter `--keep-temp`.

## Step by step guide

`quant3p` consists of four steps that can be useful by itself:
* Finding potential exons by calling peaks with MACS2,
* Fixing annotation by associating the potential exons with annotated genes,
* Fixing multimapper tags in the alignment files by selecting only alignment that maps only to annotated exons,
* Counting reads with HTSeq.

### Finding potential exons

To find potential exons we call peaks in combined RNA-seq data. We do it separetely for positive and negative strands
to use strand-specificity of our data. It's done by calling `macs2-stranded`.

```bash
macs2-stranded -n example example/bam/*.bam
```

Here:
* `-n example` sets the name of experiment to be `example`,
* `example/bam/*.bam` is a list of bam-files to process.

For the example it should produce 15 peaks (file `example_peaks.narrowPeak`):
```
chr14  25846959   25847377   example.pos_peak_1   372   +  13.33333  42.28070   37.23716   163
chr14  25871156   25871413   example.pos_peak_2   114   +  9.52381   16.27811   11.43326   80
chr14  25886465   25886948   example.pos_peak_3   1319  +  25.60976  138.30992  131.94069  273
chr14  26026235   26026718   example.pos_peak_4   1324  +  25.81967  138.85394  132.45966  273
chr14  26165849   26166332   example.pos_peak_5   1324  +  25.81967  138.85394  132.45966  273
chr5   140752996  140753373  example.pos_peak_6   718   +  22.40260  77.17718   71.88202   214
chr6   83341577   83341950   example.pos_peak_7   724   +  6.68317   77.78582   72.48792   197
chr6   83342407   83342877   example.pos_peak_8   1071  +  8.43220   112.89129  107.17482  249
chr6   83343311   83343810   example.pos_peak_9   1078  +  8.51155   113.64171  107.87556  272
chr6   83351316   83351531   example.pos_peak_10  173   +  9.23567   22.26187   17.36432   79
chr6   83353068   83353297   example.pos_peak_11  397   +  14.96815  44.79282   39.73446   80
chr6   83358160   83358528   example.pos_peak_12  1767  +  34.81308  185.13048  176.73187  183
chr5   140758332  140758782  example.neg_peak_1   1403  -  28.53982  148.70955  140.31094  182
chr5   140806999  140807270  example.neg_peak_2   73    -  7.81250   13.03361   7.39012    75
chr6   83346833   83347040   example.neg_peak_3   200   -  13.04348  25.95999   20.05134   37
```

### Fixing annotaion

Next we fix annotation by adding RNA-seq peaks as exons to the transcripts if the peak overlaps with a region of 5Kbp downstream of the transcript's 3' end.

```bash
gtf-extend -g example/mm10.slice.gtf -p example_peaks.narrowPeak -o mm10.slice.fixed.gtf --extns-out mm10.slice.extensions.gtf
```
Here:
* `-g example/mm10.slice.gtf` tells to use `example/mm10.slice.gtf ` as the annotation,
* `-p example_peaks.narrowPeak` tells to use `example_peaks.narrowPeak` as peaks,
* `-o mm10.slice.fixed.gtf` tells to put fixed annotation into `mm10.slice.fixed.gtf` file,
* `--extns-out mm10.slice.extensions.gtf` tells to put newly added exons into `mm10.slice.extensions.gtf` files.

In the example we added 4 new exons (file `mm10.slice.extensions.gtf`):
```
chr6  extension  exon  83341578   83341950   724   +  .  transcript_id  "NM_145571";  gene_id  "232157"
chr6  extension  exon  83342408   83342877   1071  +  .  transcript_id  "NM_145571";  gene_id  "232157"
chr6  extension  exon  83343312   83343810   1078  +  .  transcript_id  "NM_145571";  gene_id  "232157"
chr5  extension  exon  140758333  140758782  1403  -  .  transcript_id  "NM_010302";  gene_id  "14673"
```

### Fixing multimappers

Main idea here is that we expect RNA-seq reads to be aligned to the genes.
So, we discard alignments of multimappers to the unannotated regions and 
fix multimapping annotation (`NH` SAM field) for such reads. Some of them
can be uniquely mapped to one gene.

```bash
fix-mm -g example/mm10.slice.gtf example/bam/2h_rep1.bam -o 2h_rep1.fixed.bam
```
Here:
* `-g example/mm10.slice.gtf ` tells to use `example/mm10.slice.gtf ` as the annotation,*
* `example/bam/2h_rep1.bam` tells to fix `example/bam/2h_rep1.bam` file,
* `-o 2h_rep1.fixed.bam` tells to put fixed bam-file into `2h_rep1.fixed.bam`,

In the `2h_rep1.bam` file all the 301 multimappers can be uniquely mapped to a single position covered by an exon.

### Counting with HTSeq

Counting with HTSeq is rather straightworward:
```bash
samtools view 2h_rep1.fixed.bam | htseq-count -s yes -t exon - mm10.slice.fixed.gtf
```

Here we use the fixed bam file (`2h_rep1.fixed.bam`) and annotaion (`mm10.slice.fixed.gtf`).

Output should be like this:
```
100039246               0
11744                   112
14673                   105
211623                  0
231842                  47
232157                  1799
66039                   1
78653                   147
NM_001270503_dup1       0
NM_001270503_dup2       0
NM_207229_dup1          0
NM_207229_dup2          0
__no_feature            463
__ambiguous             0
__too_low_aQual         0
__not_aligned           0
__alignment_not_unique  289
```

If we run `htseq-count` for the original data:
```bash
samtools view -h example/bam/2h_rep1.bam | htseq-count -s yes -t exon - example/mm10.slice.gtf
```

A lot of reads will not be counted:
```
100039246               0
11744                   0
14673                   0
211623                  0
231842                  44
232157                  17
66039                   0
78653                   56
NM_001270503_dup1       0
NM_001270503_dup2       0
NM_207229_dup1          0
NM_207229_dup2          0
__no_feature            2020
__ambiguous             0
__too_low_aQual         0
__not_aligned           0
__alignment_not_unique  826
```

## Dependencies

`macs2-stranded` depends on `samtools` and `macs2` executables.  `macs2` version should be >= 2.0.10

`gtf-extend` and `fix-mm` needs `pybedtools` and `pysam` python packages and `bedtools` with `samtools` installed


All the dependencies can be installed using conda: `conda install -c bioconda macs2 pybedtools pysam HTseq`
