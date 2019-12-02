#!/usr/bin/env python

import sys
import pysam
import HTSeq
import argparse
from functools import reduce
from copy import copy
import csv

def update_if(d, f, key, value):
    if not key in d:
        d[key] = value
        return None
    old_value = d[key]
    should_update = f(old_value, value)
    if not should_update:
        return value
    d[key] = value
    return old_value

def main():
    argparser = make_argparser()
    args = argparser.parse_args()

    extension_3p = args.extension_3p

    print("Loading peaks...")
    peaks = HTSeq.GenomicArrayOfSets("auto", stranded=True)

    print("Reading annotation...")

    for (chrom, start, end, name, score, strand, fold_enrichment, pval, qval, summit) in \
            csv.reader(open(args.peaks_file), delimiter="\t"):
        iv = HTSeq.GenomicInterval(chrom, int(start), int(end), strand)
        #peak = HTSeq.GenomicFeature(name, "peak", iv)
        #peak.score = score

        peaks[iv] += (iv, name)


    def is_upstream(e1, e2):
        iv1 = e1.iv
        iv2 = e2.iv
        assert iv1.chrom == iv2.chrom
        assert iv1.strand == iv2.strand


        if iv1.strand == "+":
            return iv1.end < iv2.end
        else:
            return iv1.start > iv2.start

    last_exons = {}

    intragenic_peaks = set()
    transcript_ivs = {}

    def update_last_exon(exon):
        key = (exon.attr["transcript_id"], exon.iv.chrom, exon.iv.strand)
        return update_if(last_exons, is_upstream, key, exon)

    def update_transcript_iv(exon):
        key = (exon.attr["transcript_id"], exon.iv.chrom, exon.iv.strand)
        if not key in transcript_ivs:
            transcript_ivs[key] = exon.iv.copy()
        else:
            transcript_ivs[key].extend_to_include(exon.iv)



    gtf_out = open(args.output_file, "w")


    for feature in HTSeq.GFF_Reader(args.gtf_file):
        if feature.type == "exon":
            update_last_exon(feature)
            update_transcript_iv(feature)

        gtf_out.write(feature.get_gff_line())

    for transcript_iv in transcript_ivs.values():
        steps = peaks[transcript_iv].steps()
        values = [value for (iv, value) in steps]
        intragenic_peaks |= reduce(set.union, values, set())

    
    print("Number of intragenic peaks:", len(intragenic_peaks))
    exons_added = 0

    extns_out = None

    if args.extns_out:
        extns_out = open(args.extns_out, "w")

    for exon in last_exons.values():
        iv = exon.iv
        post3p_iv = iv.copy()
        if iv.strand == "+":
            post3p_iv.start = iv.end
            post3p_iv.end = iv.end + extension_3p
        elif iv.strand == "-":
            post3p_iv.end = iv.start
            post3p_iv.start = max(0, iv.start - extension_3p)

        if post3p_iv.length <= 0:
            continue
    
        steps = peaks[post3p_iv].steps()
        values = [value for (iv, value) in steps]
        overlapping_peaks = reduce(set.union, values, set())

        overlapping_peaks.difference_update(intragenic_peaks)

        for (iv, name) in overlapping_peaks:
            extension = HTSeq.GenomicFeature(name, "exon", iv)
            extension.source = args.source
            extension.attr = exon.attr
            exons_added += 1
            gtf_out.write(extension.get_gff_line())
            if extns_out:
                extns_out.write(extension.get_gff_line())

    gtf_out.close()
    print("Exons added:", exons_added)


def make_argparser():
    parser = argparse.ArgumentParser(description="Extend annotation with MACS peaks")

    parser.add_argument(
            '--ext-3p',
            dest="extension_3p",
            metavar="EXTENSION",
            type=int,
            default=5000,
            help="how far 3' exon end can be potentially extended (default = 5000)")

    parser.add_argument(
            '-g', '--annotation',
            dest="gtf_file",
            metavar="FILE",
            required=True,
            help="gtf/gff/... file with genome annotation")

    parser.add_argument(
            '-p', '--peaks',
            dest="peaks_file",
            metavar="FILE",
            required=True,
            help=".narrowPeak file with MACS2 peaks")

    parser.add_argument(
            '-o', '--out',
            dest="output_file",
            metavar="GTF",
            required=True,
            help="file to write extended annotation")

    parser.add_argument(
            '--extns-out',
            dest="extns_out",
            metavar="GTF",
            help="gtf file to print added exons")

    parser.add_argument(
            '--source',
            dest="source",
            metavar="NAME",
            default="extension",
            help="value to put to in the GTF source column (default = extension)")

    return parser



if __name__ == '__main__':
    main()
