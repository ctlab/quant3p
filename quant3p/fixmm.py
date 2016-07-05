#!/usr/bin/env python

import os
import sys
import pybedtools
import argparse
import pysam
import logging
from logging import debug, info, warn, error
from copy import copy


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

def is_upstream(e1, e2):
    iv1 = e1
    iv2 = e2
    assert iv1.chrom == iv2.chrom
    assert iv1.strand == iv2.strand


    if iv1.strand == "+":
        return iv1.end < iv2.end
    else:
        return iv1.start > iv2.start

def extended_exons(features, extension_5p, extension_3p, extend_all_exons):
    def extend_iv(iv):
        res = copy(iv)
        if res.strand == "+":
            res.end += extension_3p
            res.start = max(res.start - extension_5p, 0)
        elif res.strand == "-":
            res.end += extension_5p
            res.start = max(res.start - extension_3p, 0)
        return res

    last_exons = {}

    def update_last_exon(exon):
        key = (exon.attrs["transcript_id"], exon.chrom, exon.strand)
        return update_if(last_exons, is_upstream, key, exon)

    for feature in features:
        if feature.fields[2] == "exon":
            iv = feature

            if extend_all_exons:
                yield extend_iv(iv)
            else:
                feature_to_push = update_last_exon(feature)
                if feature_to_push:
                    yield feature_to_push

    for exon in last_exons.itervalues():
        yield extend_iv(exon)

def main():
    logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)-15s - %(levelname)s: %(message)s"
            )
    argparser = make_argparser()
    args = argparser.parse_args()

    extension_3p = args.extension_3p
    extension_5p = args.extension_5p

    if not args.stats_only and not args.output_file:
        error("Please either provide output file or --stats-only flag")
        sys.exit(1)

    info("Finding exonic multimappers...")

    xgtf = pybedtools.BedTool(extended_exons(
            pybedtools.BedTool(args.gtf_file),
            extension_5p=extension_5p,
            extension_3p=extension_3p,
            extend_all_exons=args.extend_all_exons))


    bamTool = pybedtools.BedTool(args.input_file)

    tocheck = {}

    exonicBam = bamTool.intersect(xgtf, split=True, s=True, nonamecheck=True)

    try:
        # using HTSeq here because of genomic intervals
        for alignment in exonicBam:
            is_multimapper = not "NH:i:1" in alignment.fields[9:]
            if not is_multimapper:
                continue

            qname = alignment.fields[0]
            if not qname in tocheck:
                tocheck[qname] = 0
            tocheck[qname] += 1
    finally:
        os.unlink(exonicBam.fn)

    info("Determining which multimappers are not real...")
    tofix = set()
    for (qname, count) in tocheck.iteritems():
        if count == 1:
            tofix.add(qname)

    if args.counts_out:
        with open(args.counts_out, "w") as counts_out:
            for (qname, count) in tocheck.iteritems():
                counts_out.write("%s\t%s\n" % (qname, count))

    info("Number of exonic mutlimappers: %d", len(tocheck))
    info("Number of fixable mutlimappers: %d", len(tofix))

    if args.stats_only:
        return

    in_samfile = pysam.Samfile(args.input_file, "rb")
    out_samfile = pysam.Samfile(args.output_file, "wb", template=in_samfile)

    for read in in_samfile:
        if read.qname in tofix:
            read.tags = [("NH", 1), ("HI", 1)] + \
                    [tag for tag in read.tags if not tag[0] in ["NH", "HI"]]
            read.mapq = args.new_aqual
        out_samfile.write(read)

    in_samfile.close()
    out_samfile.close()


def make_argparser():
    parser = argparse.ArgumentParser(description="Fix NH/HI flags for potentially exonic multimappers")

    parser.add_argument(
            '--ext-3p',
            dest="extension_3p",
            metavar="EXTENSION",
            type=int,
            default=10000,
            help="how far extend 3' exon ends (default = 10000)")

    parser.add_argument(
            '--ext-5p',
            dest="extension_5p",
            metavar="EXTENSION",
            default=0,
            type=int,
            help="how far extend 5' exon ends (defaults = 0)")

    parser.add_argument(
            '--aqual',
            dest="new_aqual",
            metavar="QUALITY",
            default=30,
            type=int,
            help="new alignment quality to set for fixed reads (default = 30)")

    parser.add_argument(
            '--all-exons',
            dest="extend_all_exons",
            default=False,
            action='store_true',
            help="Extend all exons of a transcript, not just the last one")

    parser.add_argument(
            '--only-last',
            dest="extend_all_exons",
            default=False,
            action='store_false',
            help="Extend only the last exon of a transcript (default)")

    parser.add_argument(
            '-g', '--annotation',
            dest="gtf_file",
            metavar="FILE",
            required=True,
            help="gtf/gff/... file with genome annotation")

    parser.add_argument(
            '--stats-only',
            dest="stats_only",
            default=False,
            action='store_true',
            help="Just show multimappers statistics, do not write the output bam")

    parser.add_argument(
            '--counts-out',
            dest="counts_out",
            metavar="FILE",
            help="file to print multimappers exonic coutns")

    parser.add_argument(
            '-o', '--out',
            dest="output_file",
            metavar="BAM",
            help="file to write fixed alignments")

    parser.add_argument(
            'input_file',
            metavar="INPUT_BAM",
            help="file with alignments to fix")
    return parser


if __name__ == '__main__':
    main()
