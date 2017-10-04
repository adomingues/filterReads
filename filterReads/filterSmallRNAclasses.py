#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import (absolute_import, division,
                        print_function)
import argparse
import pysam
import sys

usage = '''
   Filter classes of small RNAs from library.
   Classes are defined by read sequence length
   '''


def getArgs():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser(
        description=usage,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '-i', '--inputBam',
        required=True,
        type=str,
        help='Path to alignments (bam) to be filtered. The bam file needs to be indexed.'
    )

    parser.add_argument(
        '-o', '--outputBam',
        required=True,
        type=str,
        help='Output Bam file to store filtered results. Use "stdout" to write output to standard out'
    )

    parser.add_argument(
        '-c', '--RNAclass',
        required=False,
        type=str,
        help='Pre-defined class of small RNAs to be output. Options: 21U, 22G or 26G'
    )

    parser.add_argument(
        '-l', '--length',
        required=False,
        type=int,
        help='The exact length (l) of reads to be kept. Not implemented'
    )

    parser.add_argument(
        '-m', '--min',
        required=False,
        type=int,
        default=0,
        help='Reads with length >= m will be keep. Not implemented'
    )

    parser.add_argument(
        '-M', '--max',
        required=False,
        type=int,
        default=100000,
        help='Reads with length <= M will be keep. Not implemented'
    )

    parser.add_argument(
        '-n', '--nuc',
        required=False,
        type=int,
        help='Alignments with n at the first position, or the reverse complement at the the last position if read is mapped in the reverse strand, will be kept. Not implemented'
    )

    args = parser.parse_args()
    return args


def eprint(*args, **kwargs):
    """
    helper function to print to stderr
    source:
    https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
    """
    print(*args, file=sys.stderr, **kwargs)


if __name__ == '__main__':
    args = getArgs()
    in_file = args.inputBam
    out_file = args.outputBam

    if args.RNAclass:
        eprint("Using preset settings for small RNA class: %s"
              % (args.RNAclass))
        if args.RNAclass == "21U":
            length = 21
            nucleotide_for = "T"
            nucleotide_rev = "A"
        elif args.RNAclass == "22G":
            min_length = 22
            max_length = 23
            nucleotide_for = "G"
            nucleotide_rev = "C"
        elif args.RNAclass == "26G":
            length = 26
            nucleotide_for = "G"
            nucleotide_rev = "C"


    piRNA_reads = 0
    other_reads = 0

    inbam = pysam.AlignmentFile(in_file, "rb")
    if out_file == "stdout":
        outbam = pysam.AlignmentFile("-", "wb", template=inbam)
    else:
        outbam = pysam.AlignmentFile(out_file, 'wb', template=inbam)

    for read in inbam.fetch():
        if (read.is_reverse is True and
                read.seq.endswith(nucleotide_rev) and
                min_length <= read.qlen <= max_length):
            piRNA_reads = piRNA_reads + 1
            outbam.write(read)
        if (read.is_reverse is False and
                read.seq.startswith(nucleotide_for) and
                min_length <= read.qlen <= max_length):
            piRNA_reads = piRNA_reads + 1
            outbam.write(read)
        else:
            other_reads = other_reads + 1

    inbam.close()

    if out_file != "stdout":
        outbam.close()
        pysam.index(out_file)

    eprint("Number of %s reads: %d" % (args.RNAclass, piRNA_reads))
    eprint("Number of reads filtered out: %d" % other_reads)
