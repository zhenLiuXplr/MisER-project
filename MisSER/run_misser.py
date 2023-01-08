#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import textwrap as _textwrap
from .fix_small_exon_multiproc import run_multiprocess
from .version import VERSION

class CustomHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string

    def _fill_text(self, text, width, indent):
        return ''.join(indent + line for line in text.splitlines(keepends=True))

    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        return _textwrap.wrap(text, 80)


"""
The script is used to import MisSER package. 
Run run_multiprocess function to find and fix small-exon-missed reads.
"""


def main():
    parser = argparse.ArgumentParser(
        description="Find and fix small-exon-missed region in noisy-long reads based on transcript annotation.",
        formatter_class=CustomHelpFormatter
    )
    group = parser.add_mutually_exclusive_group()
    # positional arguments
    parser.add_argument("inBam", help="Input original bam file.")
    parser.add_argument(
        "genomeFasta", 
        help="Reference genome fasta file, with fai index under the same directory."
        )
    parser.add_argument("annotBed", help="Annotated transcripts file in BED12 format.")
    parser.add_argument("outRegion", help="Output Region file, regions contain missed small exons.")
    # optional arguments
    parser.add_argument(
        "-v", "--version", help="Print version and exit.",
        action="version", version="MisSER {0}".format(VERSION)
    )
    parser.add_argument(
        "-c", "--coreNum", help="The number of cpu cores we used.",
        default=1, type=int, metavar="N"
    )
    parser.add_argument(
        "-s", "--exonSizeThd", 
        help="The threshold of exons size, ignore exons with size > exonSizeThd.",
        default=80, type=int, metavar="N"
    )
    parser.add_argument(
        "-d", "--deltaRatioThd",
        help="The threshold of absolute delta ratio, ignore abs(delta ratio) > deltaRatioThd.",
        default=0.5, type=float, metavar="float"
    )
    parser.add_argument(
        "-f", "--flankLen", help="The extended length on the both sides of realign region.",
        default=20, type=int, metavar="N"
    )
    parser.add_argument(
        "--strandSpecific", help="Only compare reads consistent with annotated strand.",
        action="store_true", default=False
    )
    parser.add_argument(
        "--allTranscripts", help="Return all possible missed exons among overlapping transcripts.",
        action="store_true", default=False
    )
    parser.add_argument(
        "--fixFlankLen", help="Flank length will not extend to cover the adjacent indels.",
        action="store_true", default=False
    )
    parser.add_argument(
        "--debugMode", help="Won't stop when meet an error in one read.",
        action="store_true", default=False
    )
    parser.add_argument(
        "--setTag", help="Set fr tags on the realigned reads",
        action="store_true", default=False
    )
    group.add_argument("-o", "--outBam", help="Output realigned bam file.", metavar="file")
    group.add_argument(
        "--onlyRegion", help="Only return the Region file without realign process.",
        action="store_true", default=False
    )
    # output help message when no param.
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    if not args.onlyRegion:
        if args.outBam is None:
            raise Exception(
                "Require -o or --outBam if --onlyRegion flag is not set!")
    run_multiprocess(
        args.inBam, args.genomeFasta, args.annotBed, args.outBam, args.outRegion,
        small_exon_size=args.exonSizeThd, flank_len=args.flankLen,
        ignore_strand=not(args.strandSpecific), nprocess=args.coreNum,
        delta_ratio_thd=args.deltaRatioThd, simplify=not(args.allTranscripts),
        float_flank_len=not(args.fixFlankLen), only_region=args.onlyRegion,
        debug_mode=args.debugMode, set_tag=args.setTag
    )
