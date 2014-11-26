#!/usr/bin/env python

import sys
import HTSeq
import argparse

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

def main(args):
    exons = HTSeq.GenomicArray("auto", stranded=True, typecode='b')

    extension_3p = args.extension_3p
    extension_5p = args.extension_5p

    if not args.stats_only and not args.output_file:
        print "Please either proved output file or --stats-only flag"
        sys.exit(1)

    print "Reading annotation..."

    def extend_iv(iv):
        res = iv.copy()
        if res.strand == "+":
            res.end += extension_3p
            res.start -= extension_5p
        elif res.strand == "-":
            res.end += extension_5p
            res.start -= extension_3p
        res.start = max(res.start, 0)
        return res

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

    def update_last_exon(exon):
        key = (exon.attr["transcript_id"], exon.iv.chrom, exon.iv.strand)
        return update_if(last_exons, is_upstream, key, exon)


    for feature in HTSeq.GFF_Reader(args.gtf_file):
        if feature.type == "exon":
            iv = feature.iv

            if args.extend_all_exons:
                iv = extend_iv(iv)
            else:
                feature_to_push = update_last_exon(feature)
                if feature_to_push:
                    iv = feature_to_push.iv
                else:
                    iv = None

            if iv and iv.length > 0:
                exons[iv] = True


    for exon in last_exons.itervalues():
        iv = exon.iv
        iv = extend_iv(iv)
        if iv.length > 0:
            exons[iv] = True


    print "Finding exonic multimappers..."
    tocheck = {}

    # using HTSeq here because of genomic intervals

    for alignment in HTSeq.BAM_Reader(args.input_file):
        is_multimapper = alignment.optional_field( "NH" ) > 1
        if not is_multimapper:
            continue

        should_check = False

        iv_seq = ( co.ref_iv for co in alignment.cigar if co.type == "M" and co.size > 0 )

        for genomic_iv in iv_seq:
            for iv, value in exons[genomic_iv].steps():
                if value:
                    should_check = True
                    break
            if should_check:
                break

        if should_check:
            qname = alignment.read.name
            if not qname in tocheck:
                tocheck[qname] = 0
            tocheck[qname] += 1

    print "Determining which multimappers are not real..."
    tofix = set()
    for (qname, count) in tocheck.iteritems():
        if count == 1:
            tofix.add(qname)

    if args.counts_out:
        with open(args.counts_out, "w") as counts_out:
            for (qname, count) in tocheck.iteritems():
                counts_out.write("%s\t%s\n" % (qname, count))

    print "Number of exonic mutlimappers:", len(tocheck)
    print "Number of fixable mutlimappers:", len(tofix)

    if args.stats_only:
        return

    infile = HTSeq.BAM_Reader(args.input_file)
    outfile = HTSeq.BAM_Writer.from_BAM_Reader(args.output_file, infile)

    for aln in infile:
        if aln.read.name in tofix:
            aln.optional_fields = [("NH", 1), ("HI", 1)] + \
                    [tag for tag in aln.optional_fields if tag[0] in ["NH", "HI"]]
        outfile.write(aln)

    outfile.close()


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
    argparser = make_argparser()
    args = argparser.parse_args()
    main(args)
