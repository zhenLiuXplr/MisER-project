#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BAM file will be input use pysam and each line in BAM file will be looped to find the possible missed small exons, with
    the threshold: small_exon_size <= 100 (default) and abs(delta_ratio) <= 0.5 (default).
The realign region will only be aligned to the exon sequences of the transcript, not compare with intron sequences.
"""

import csv
import re
import sys
import time
import warnings
import os
import logging

import numpy as np
import pandas as pd
import multiprocessing as mp
import pysam
import parasail
import collections

from itertools import groupby
from .block_class import (Cigar, Cigar_fromlist, ReadBlock, CigarSliceError)
from .find_small_exon import find_missed_exons


def groupby_unsorted(seq, key_fun=lambda x: x):
    """
    Unsorted version of groupby in itertools
    """
    indexs_dict = collections.defaultdict(list)
    for i, elem in enumerate(seq):
        indexs_dict[key_fun(elem)].append(i)
    for key, idxs in indexs_dict.items():
        yield key, [seq[i] for i in idxs]


def parse_bed(bed_path):
    """
    The function parse an annotation BED file to a dataframe contain list.
    We used to use list of dict as the result.
    [{'tx_i':146, 'chrom':'chrIS', 'strand':'+', 'start':6955486, 'end':6962816, 'name':'R1_101_1', 'block_count':3, 'block_sizes':[243, 116, 360], 'block_starts':[0, 4897, 6970]},
     {...},
     ...
     {...}]
    Actually, the list of dict is easy to be transformed to dataframe by pd.DataFrame(list_dict) and vice versa df.to_dict('records')
    The result dataframe should be sorted by start and end sites.
    """
    bed_col = [
        "chrom", "chromStart", "chromEnd", "name", "score", "strand",
        "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"
    ]
    # Set the chrom as str because reference in ensembl do not have prefix 'chr'.
    # the read.reference_name in pysam is string.
    bed = pd.read_csv(
        bed_path, sep='\t', header=None, names=bed_col,
        quoting=csv.QUOTE_NONE,
        dtype={"chrom": str, "blockSizes": str, "blockStarts": str}
    )
    bed["tx_i"] = bed.index
    #
    bed_dict_list = bed.to_dict("records")
    parse_dict_list = []
    for tx_dict in bed_dict_list:
        tx_dict["blockSizes"] = [np.int64(x) for x in tx_dict["blockSizes"].split(',') if len(x) > 0]
        tx_dict["blockStarts"] = [np.int64(x) for x in tx_dict["blockStarts"].split(',')if len(x) > 0]
        tx_dict["blockNum"] = len(tx_dict["blockSizes"])
        parse_dict_list.append(tx_dict)
    bed = pd.DataFrame(parse_dict_list)
    bed_df = pd.DataFrame(bed[["tx_i", "chrom", "strand"]])
    bed_df["tx"] = bed["name"]
    bed_df["start"] = bed["chromStart"]
    bed_df["end"] = bed["chromEnd"]
    bed_df["exon_sizes"] = bed["blockSizes"]
    bed_df["exon_starts"] = bed["blockStarts"]
    bed_df["exon_num"] = bed["blockNum"]
    return bed_df


def extend_tx(tx):
    """
    The function extends the list at exon_sizes and exon_starts in the transcript dataframe into exons.
    One line stands for one exon in the extended result.
    """
    if not all(tx["start"] < tx["end"]):  # We require start < end in the bed_df
        raise Exception("Error: start ≥ end in bed file!")
    exon_starts = tx.apply(
        lambda r: r["start"] + pd.Series(r["exon_starts"]), axis=1)
    exon_ends = tx.apply(
        lambda r: r["start"] + pd.Series(r["exon_starts"]) + pd.Series(r["exon_sizes"]), axis=1)
    exon_starts = exon_starts.stack().reset_index(level=1)
    exon_ends = exon_ends.stack().reset_index(level=1)
    if (not all(exon_starts["level_1"] == exon_ends["level_1"])) or (not all(exon_starts.index == exon_ends.index)):
        raise Exception("Error: extend exon fall!")
    exons = pd.DataFrame(
        {"exon_i": exon_starts["level_1"], "exon_start": exon_starts.iloc[:, 1], "exon_end": exon_ends.iloc[:, 1]}, dtype="int")
    exons = tx.drop(["exon_sizes", "exon_starts"], axis=1).join(exons)
    exons.set_index(["tx_i", "exon_i"], inplace=True)
    return exons


def trace_realign_region(row, read_cigar, read_start):
    """
    The function is used to trace the rdp and ci at the border of realign region and return a dict
    """
    realn_start_t = read_cigar.trace_ci_rdp(start_ci=(0, 0), start_pos=read_start,
                                            end_pos=row["realign_start"], start_rdp=0)
    if not realn_start_t[2]:
        raise Exception("Error: cannot trace to the realign border", row)
    realn_end_t = read_cigar.trace_ci_rdp(start_ci=realn_start_t[0], start_pos=row["realign_start"],
                                          end_pos=row["realign_end"], start_rdp=realn_start_t[1])
    if not realn_end_t[2]:
        raise Exception("Error: cannot trace to the realign border", row)
    return {"realign_start_ci": realn_start_t[0], "realign_end_ci": realn_end_t[0],
            "realign_start_rdp": realn_start_t[1], "realign_end_rdp": realn_end_t[1]}


def edit_margin_len(row, read_cigar):
    """
    The function is used to edit the margin length in the original realign_ss.
    Margin length is called margin_length in the old file.
    One insertion 'I' or 'D' will +1 or -1 in the margin length.
    """
    cigar_slice = read_cigar.cslice(
        start_ci=row["realign_start_ci"], end_ci=row["realign_end_ci"])
    new_margin_len = row["margin_len"]
    for n, c in cigar_slice:
        if c == 'D':
            new_margin_len -= n
        elif c == 'I':
            new_margin_len += n
    return new_margin_len


def extend_realign_region(row, read_cigar):
    """
    This part is separated from the original modify_realign_region.

    Because minimap use bonus strategy on the splice site, There might be a long deletion or insertion beside the
    splice site. In case to make sure that the realign start do not begin with a deletion or the realign end do not
    end with a insertion and deletion, we need to extend the boundary of the realign region when the case happen.

    Shown as follows:

    Deletion:
    Start:  ... End:
      # @         @ #
    ATGCA       GCTTC
    AT--A       GC--C
        @           @

    Insertion:
    Start:  ... End:
        @           @
    AT--A       GC--C
    ATGCA       GCTTC
      # @         @ #

    The symbol '@' means the original position of ref_start/end (above the sequence) and query_start/end (below the
    sequence); The symbol '#' means the position after modification.
    """
    try:
        left_slice_group = read_cigar.cslice(
            end_ci=row["realign_start_ci"])[-1]
        if left_slice_group[1] in {'D', 'I'}:
            if left_slice_group[1] == 'D':
                row["realign_start"] -= left_slice_group[0]
            else:
                row["realign_start_rdp"] -= left_slice_group[0]
            if row["realign_start_ci"][1] == 0:
                row["realign_start_ci"] = (row["realign_start_ci"][0]-1, 0)
            else:
                row["realign_start_ci"] = (row["realign_start_ci"][0], 0)
    except IndexError:  # most case because realign_start_ci is (0, 0)
        pass
    try:
        right_slice_group = read_cigar.cslice(
            start_ci=row["realign_end_ci"])[0]
        if right_slice_group[1] in {'D', 'I'}:
            if right_slice_group[1] == 'D':
                row["realign_end"] += right_slice_group[0]
            else:
                row["realign_end_rdp"] += right_slice_group[0]
            row["realign_end_ci"] = (row["realign_end_ci"][0]+1, 0)
    except IndexError:  # most case because realign_end_ci reach the cigar end.
        pass
    return row


def exclude_clip(row, read_cigar):
    """
    Separated from extend_realign_region, used to exclude the soft and hard clip at the start or end of the query_sequence.
    Because sometimes realign region might include the softclip region and the softclip region usually is really
    long which can lead to a wrong realignment result.
    Realignment will overwrite the hard clip, so we need to exclude them.
    """
    # Exclude the 'S' at both ends of the realign region
    cigar_slice = read_cigar.cslice(
        start_ci=row["realign_start_ci"], end_ci=row["realign_end_ci"])
    if cigar_slice[0][1] == 'S':
        row["realign_start_rdp"] += cigar_slice[0][0]
        row["realign_start_ci"] = (row["realign_start_ci"][0]+1, 0)
    if cigar_slice[0][1] == 'H':
        row["realign_start_ci"] = (row["realign_start_ci"][0]+1, 0)
    if cigar_slice[-1][1] == 'S':
        row["realign_end_rdp"] -= cigar_slice[0][0]
        row["realign_end_ci"] = (row["realign_end_ci"][0], 0)
    if cigar_slice[-1][1] == 'H':
        row["realign_end_ci"] = (row["realign_end_ci"][0], 0)
    return row


def modify_realign_region(row, read_cigar, read_start, float_flank_len):
    """
    The function is used to modify the realignment region.
    It also combines the trace, edit margin length and extend realign region process.
    Return: modified realign_ss
    """
    # Tracing process
    trace_dict = trace_realign_region(row[["realign_start", "realign_end", "read_intron_start",
                                           "read_intron_end"]], read_cigar, read_start)
    row = row.append(pd.Series(trace_dict))
    if float_flank_len:
        # Extend the boundary of realign region
        row = extend_realign_region(row, read_cigar)
    # Exclude soft clip
    row = exclude_clip(row, read_cigar)
    # Editing margin length
    row["margin_len_mod"] = edit_margin_len(row, read_cigar)
    row["delta_ratio_mod"] = (
        row["margin_len_mod"] - row["sum_exon_size"]) / row["sum_exon_size"]
    return row


def trace_realign_result(row, realign_cigar):
    """
    The function is used to trace the cigar index of introns which are used to split the realign_cigar.
    Return a ci_list.
    realn_tx_pos means the position regardless of introns and relative to the realign_start. The position of
    bases on the tx_seq.
    """
    ci_list = []
    intron_start = row["tx_leftexon_end"]
    realn_tx_pos = row["tx_leftexon_end"]-row["realign_start"]
    realn_intron_pos = [realn_tx_pos]
    for exon_size in row["exon_sizes"]:
        realn_tx_pos += exon_size
        realn_intron_pos.append(realn_tx_pos)
    for realn_tx_pos in realn_intron_pos:
        trace = realign_cigar.trace_ci_rdp(
            start_ci=(0, 0), start_pos=0, end_pos=realn_tx_pos, start_rdp=0)
        if not trace[2]:
            raise Exception(
                "Error:, cannot trace to the intron position in the realign cigar!")
        ci_list.append(trace[0])
    return ci_list


def cal_alignment_score(cigar_list, read_seq, ref_seq, cm, mm, go, ge):
    """
    mm, go and ge usually are negative to give penalty.
    """
    score = 0
    pos = 0
    rdp = 0
    for n, c in cigar_list:
        if n < 1:
            raise Exception(
                "Error: list in cigar is not standard!", cigar_list)
        if c == 'M':
            for _ in range(n):
                if read_seq[rdp] == ref_seq[pos]:
                    score += cm
                else:
                    score += mm
                pos += 1
                rdp += 1
        elif c == '=':
            score += cm*n
            pos += n
            rdp += n
        elif c == 'X':
            score += mm*n
            pos += n
            rdp += n
        elif c in {'I'}:
            score += go + ge*(n-1)
            rdp += n
        elif c == 'D':
            score += go + ge*(n-1)
            pos += n
        elif c == 'N':
            pos += n  # we do not give penalty on introns
        elif c == 'P':
            warnings.warn("Warning: 'P' in the cigar!", c)
        else:
            raise Exception(
                "Error: cigar character exceed assumption {'M', '=', 'X', 'I', 'D', 'P', 'N'}!", c)
    if rdp != len(read_seq) or pos != len(ref_seq):
        raise Exception("Error: input cigar not correspond to read_seq or ref_seq!", '\n',
                        cigar_list, '\n',
                        read_seq, '\n',
                        ref_seq)
    return score


def check_cigar_seq(cigar, read_seq, ref_seq):
    """
    The function is used to check whether the cigar is corresponded to read_seq and ref_seq
    """
    pos = 0
    rdp = 0
    for n, c in cigar.list:
        if c in {'M', 'X', '='}:
            pos += n
            rdp += n
        elif c in {'S', 'I'}:
            rdp += n
        elif c in {'D', 'N'}:
            pos += n
    if rdp != len(read_seq) or pos != len(ref_seq):
        raise Exception("Error: input cigar not correspond to read_seq or ref_seq!", '\n',
                        cigar.list, '\n',
                        read_seq, '\n',
                        ref_seq)
    return 0


def realign_read(read, read_cigar, realign_ss, ref_fa, score_matrix):
    """
    The function realign the region indicated in the realign_ss.
    Suffix rel means relative. pos_rel means the reference position relative to the realign start site.
    pos_rel = pos - realign_start.
    """
    check_cigar_seq(read_cigar, read.query_sequence, ref_fa.fetch(
        realign_ss["chrom"], start=read.reference_start, end=read.reference_end))
    realign_start = realign_ss["realign_start"]
    realign_end = realign_ss["realign_end"]
    # fetch half-open interval.
    # Use upper() to avoid soft-masked sequence with lowercase
    ref_seq = ref_fa.fetch(
        realign_ss["chrom"], start=realign_start, end=realign_end).upper()
    # tx_seq: part of ref_seq which does not contain intron sequences.
    # Also calculate transcript intron lengths
    intron_lens = []
    intron_start = realign_ss["tx_leftexon_end"]
    tx_seq = ref_seq[:realign_ss["tx_leftexon_end"]-realign_start]
    for s, e in zip(realign_ss["small_exon_starts"], realign_ss["small_exon_ends"]):
        tx_seq += ref_seq[s-realign_start:e-realign_start]
        intron_lens.append(s - intron_start)
        intron_start = e
    intron_lens.append(realign_ss["tx_rightexon_start"] - intron_start)
    tx_seq += ref_seq[realign_ss["tx_rightexon_start"]-realign_start:]
    read_seq = read.query_sequence[realign_ss["realign_start_rdp"]:realign_ss["realign_end_rdp"]]
    realign_result = parasail.nw_trace(read_seq, tx_seq, 1, 1, score_matrix)
    realign_score = realign_result.score
    realign_cigar_str = realign_result.cigar.decode.decode("utf-8")
    # We find the parasail might confuse the "X" with "=", so we convert all "X" and "=" to "M"
    realign_cigar_str = re.sub("[X=]", "M", realign_cigar_str)
    realign_cigar = Cigar(realign_cigar_str)
    realign_cigar.standardize()
    # check cal_alignment_score
    if realign_score != cal_alignment_score(realign_cigar.list, read_seq, tx_seq, 1, -1, -1, -1):
        print("read_seq:{0}\ntx_seq:{1}\nrealign_cigar_str:{2}\n".format(
            read_seq, tx_seq, realign_cigar_str
        ))
        raise Exception("Error: parasail alignment score different with the calculation!", realign_score,
                        cal_alignment_score(realign_cigar.list, read_seq, tx_seq, 1, -1, -1, -1))
    ori_realn_cigar_list = read_cigar.cslice(start_ci=realign_ss["realign_start_ci"],
                                             end_ci=realign_ss["realign_end_ci"])
    ori_score = cal_alignment_score(
        ori_realn_cigar_list, read_seq, ref_seq, 1, -1, -1, -1)
    increase_score = realign_score - ori_score
    realn_intron_cis = trace_realign_result(
        realign_ss[["tx_leftexon_end", "realign_start", "exon_sizes"]], realign_cigar)
    # joint new cigar
    new_cigar_list = read_cigar.cslice(end_ci=realign_ss["realign_start_ci"])
    start_ci = (0, 0)
    for ci, n in zip(realn_intron_cis, intron_lens):
        new_cigar_list.extend(realign_cigar.cslice(
            start_ci=start_ci, end_ci=ci))
        new_cigar_list.extend([[n, 'N']])
        start_ci = ci
    new_cigar_list.extend(realign_cigar.cslice(start_ci=start_ci))
    new_cigar_list.extend(read_cigar.cslice(
        start_ci=realign_ss["realign_end_ci"]))
    new_cigar = Cigar_fromlist(new_cigar_list)
    # check new cigar:
    # check_cigar_seq(new_cigar, read.query_sequence, ref_fa.fetch(realign_ss["chrom"], start=read.reference_start, end=read.reference_end))
    return new_cigar, increase_score


def find_fix(proc_region, ori_bam_path, out_bam_path, ref_fa_path, annot_bed, annot_exon,
             small_exon_size=80, flank_len=20, ignore_strand=True, set_tag=True,
             delta_ratio_thd=0.5, simplify=True, float_flank_len=False, only_region=False):
    """
    The function is used to find and fix small exon in one process.
    The proc_region is a three elements dict. {proc_id, start_read_id, end_read_id}
    Return a dict {"filter_list":[a list of realign_ss], "stat_dict":{a dict with statistics of misaligned reads}}
    """
    proc_id = proc_region["proc_id"]
    score_matrix = parasail.matrix_create("ACGT", 1, -1)
    ori_bam = pysam.AlignmentFile(ori_bam_path, "rb")
    if not os.path.exists("{0}.fai".format(ref_fa_path)):
        print("Build reference index...")
        pysam.faidx(ref_fa_path)
    ref_fa = pysam.FastaFile(ref_fa_path)
    # 0 based, left close right open
    start_read_id = proc_region["start_read_id"]
    end_read_id = proc_region["end_read_id"]
    if not only_region:
        out_bam_path_proc = out_bam_path + ".tmp{0}".format(proc_id)
        realign_bam = pysam.AlignmentFile(
            out_bam_path_proc, "wb", template=ori_bam)
    filter_list = []
    read_i = 0  # 0 based
    spliced_read_num = 0
    total_intron = 0
    misaligned_read_num = 0
    misaligned_intron_num = 0
    fix_read_num = 0
    fix_intron_num = 0
    for read in ori_bam.fetch(until_eof=True):
        if not(start_read_id <= read_i and read_i < end_read_id):
            read_i += 1
            continue
        if read.cigarstring is None or read.query_sequence is None:
            if not only_region:
                realign_bam.write(read)
            read_i += 1
            if read.cigarstring is not None and "N" in read.cigarstring:
                spliced_read_num += 1
            continue
        read_miss_exon_flag = False
        read_fix_flag = False
        read_block = ReadBlock(read, read_i)
        realign_list = find_missed_exons(read_block, annot_bed, annot_exon, ignore_strand=ignore_strand,
                                         flank_len=flank_len, small_exon_size=small_exon_size)
        if realign_list is not False:
            read_cigar = read_block.cigar
            read_start = read_block.starts[0]
            fix_region_num = 0  # count for fr tags
            for ith_intron, realign_intron_group in groupby_unsorted(realign_list, lambda x: x["read_ith_intron"]):
                realn_scores_increase = []
                realn_cigars = []
                delta_ratio_mods_abs = []
                intron_miss_exon_flag = False
                intron_fix_flag = False
                for i in range(len(realign_intron_group)):
                    realign_intron_group[i] = modify_realign_region(
                        realign_intron_group[i], read_cigar, read_start, float_flank_len)
                    delta_ratio_mods_abs.append(
                        abs(realign_intron_group[i]["delta_ratio_mod"]))
                if simplify:
                    min_ind = delta_ratio_mods_abs.index(
                        min(delta_ratio_mods_abs))
                    realign_intron_group = [realign_intron_group[min_ind]]
                for realign_ss in realign_intron_group:
                    # filter with delta_ratio
                    if abs(realign_ss["delta_ratio_mod"]) <= delta_ratio_thd:
                        read_miss_exon_flag = True
                        intron_miss_exon_flag = True
                        realign_ss["delta_ratio_flag"] = True
                        # add delta_ratio_flag when delta_ratio_mod <= delta_ratio_thd
                        if not only_region:
                            realign_result = realign_read(
                                read, read_cigar, realign_ss, ref_fa, score_matrix)
                            realign_ss["realn_increase_score"] = realign_result[1]
                            realign_ss["realn_flag"] = (True if realign_result[1] > 0 else False)
                            # add realn_flag when the region is fixed
                            realn_cigars.append(realign_result[0])
                            realn_scores_increase.append(realign_result[1])
                        filter_list.append(realign_ss)
                if intron_miss_exon_flag:
                    misaligned_intron_num += 1
                    if not only_region:
                        # when realigned region improve the alignment score, we will keep the realignment with max score.
                        if max(realn_scores_increase) > 0:
                            read_cigar = realn_cigars[realn_scores_increase.index(
                                max(realn_scores_increase))]
                            # check cigar string
                            # read_cigar.getstr()
                            read_fix_flag = True
                            intron_fix_flag = True
                            fix_intron_num += 1
                            fix_region_num += 1
            # if not only_region:
            #     new_cigar_str = read_cigar.getstr()
            #     read.cigarstring = new_cigar_str
        if not only_region:
            if read_fix_flag:
                new_cigar_str = read_cigar.getstr()
                read.cigarstring = new_cigar_str
                if set_tag:
                    if read.has_tag('fr'):
                        warnings.warn("Overwrite the existed fr tags on read:", read.qname)
                    read.set_tag('fr', fix_region_num, 'i')
            realign_bam.write(read)
        read_i += 1
        if read.cigarstring is not None and "N" in read.cigarstring:
            spliced_read_num += 1
        total_intron += read_block.introns["intron_num"]
        if read_miss_exon_flag:
            misaligned_read_num += 1
        if read_fix_flag:
            fix_read_num += 1
    stat_dict = {"read_count": read_i, "spliced_read_num": spliced_read_num, "total_intron": total_intron,
                 "misaligned_read_num": misaligned_read_num,
                 "misaligned_intron_num": misaligned_intron_num,
                 "fix_read_num": fix_read_num, "fix_intron_num": fix_intron_num}
    # The read_count in stat_dict is the count of the whole bam file, not the part in the process.
    res_dict = {"filter_list": filter_list, "stat_dict": stat_dict}
    ori_bam.close()
    if not only_region:
        realign_bam.close()
        res_dict["out_bam_path_proc"] = out_bam_path_proc
    return res_dict


def find_fix_debug(proc_region, ori_bam_path, out_bam_path, ref_fa_path, annot_bed, annot_exon,
                   small_exon_size=80, flank_len=20, ignore_strand=True, set_tag=True,
                   delta_ratio_thd=0.5, simplify=True, float_flank_len=False, only_region=False):
    """
    The function is used to find and fix small exon in one process.
    The proc_region is a three elements dict. {proc_id, start_read_id, end_read_id}
    Return a dict {"filter_list":[a list of realign_ss], "stat_dict":{a dict with statistics of misaligned reads}}
    """
    proc_id = proc_region["proc_id"]
    score_matrix = parasail.matrix_create("ACGT", 1, -1)
    ori_bam = pysam.AlignmentFile(ori_bam_path, "rb")
    # define a bam to record the reads with error
    error_bam_path = ori_bam_path+".tmp{0}".format(proc_id)+"_err"
    error_bam = pysam.AlignmentFile(
        error_bam_path, "wb", template=ori_bam)
    error_case = 0
    # Check the exist of reference fasta index
    if not os.path.exists("{0}.fai".format(ref_fa_path)):
        print("Build reference index...")
        pysam.faidx(ref_fa_path)
    ref_fa = pysam.FastaFile(ref_fa_path)
    # 0 based, left close right open
    start_read_id = proc_region["start_read_id"]
    end_read_id = proc_region["end_read_id"]
    if not only_region:
        out_bam_path_proc = out_bam_path + ".tmp{0}".format(proc_id)
        realign_bam = pysam.AlignmentFile(
            out_bam_path_proc, "wb", template=ori_bam)
    filter_list = []
    read_i = 0  # 0 based
    spliced_read_num = 0
    total_intron = 0
    misaligned_read_num = 0
    misaligned_intron_num = 0
    fix_read_num = 0
    fix_intron_num = 0
    for read in ori_bam.fetch(until_eof=True):
        try:
            if not(start_read_id <= read_i and read_i < end_read_id):
                read_i += 1
                continue
            if read.cigarstring is None or read.query_sequence is None:
                if not only_region:
                    realign_bam.write(read)
                read_i += 1
                if read.cigarstring is not None and "N" in read.cigarstring:
                    spliced_read_num += 1
                continue
            read_miss_exon_flag = False
            read_fix_flag = False
            read_block = ReadBlock(read, read_i)
            realign_list = find_missed_exons(
                read_block, annot_bed, annot_exon, ignore_strand=ignore_strand,
                flank_len=flank_len, small_exon_size=small_exon_size
            )
            if realign_list is not False:
                read_cigar = read_block.cigar
                read_start = read_block.starts[0]
                fix_region_num = 0
                for ith_intron, realign_intron_group in groupby_unsorted(realign_list, lambda x: x["read_ith_intron"]):
                    realn_scores_increase = []
                    realn_cigars = []
                    delta_ratio_mods_abs = []
                    intron_miss_exon_flag = False
                    intron_fix_flag = False
                    for i in range(len(realign_intron_group)):
                        realign_intron_group[i] = modify_realign_region(
                            realign_intron_group[i], read_cigar, read_start, float_flank_len)
                        delta_ratio_mods_abs.append(
                            abs(realign_intron_group[i]["delta_ratio_mod"]))
                    if simplify:
                        min_ind = delta_ratio_mods_abs.index(
                            min(delta_ratio_mods_abs))
                        realign_intron_group = [realign_intron_group[min_ind]]
                    for realign_ss in realign_intron_group:
                        # filter with delta_ratio
                        if abs(realign_ss["delta_ratio_mod"]) <= delta_ratio_thd:
                            read_miss_exon_flag = True
                            intron_miss_exon_flag = True
                            realign_ss["delta_ratio_flag"] = True
                            # add delta_ratio_flag when the region is found
                            if not only_region:
                                realign_result = realign_read(
                                    read, read_cigar, realign_ss, ref_fa, score_matrix)
                                realign_ss["realn_increase_score"] = realign_result[1]
                                realign_ss["realn_flag"] = (True if realign_result[1] > 0 else False)
                                # add realn_flag when the region is fixed
                                realn_cigars.append(realign_result[0])
                                realn_scores_increase.append(realign_result[1])
                            filter_list.append(realign_ss)
                    if intron_miss_exon_flag:
                        misaligned_intron_num += 1
                        if not only_region:
                            # when realigned region improve the alignment score, we will keep the realignment with max score.
                            if max(realn_scores_increase) > 0:
                                read_cigar = realn_cigars[realn_scores_increase.index(
                                    max(realn_scores_increase))]
                                # check cigar string
                                # read_cigar.getstr()
                                read_fix_flag = True
                                intron_fix_flag = True
                                fix_intron_num += 1
                                fix_region_num += 1
                # if not only_region:
                #     new_cigar_str = read_cigar.getstr()
                #     read.cigarstring = new_cigar_str
            if not only_region:
                if read_fix_flag:
                    new_cigar_str = read_cigar.getstr()
                    read.cigarstring = new_cigar_str
                    if set_tag:
                        if read.has_tag('fr'):
                            warnings.warn("Overwrite the existed fr tags on read:", read.qname)
                        read.set_tag('fr', fix_region_num, 'i')
                realign_bam.write(read)
            read_i += 1
            if read.cigarstring is not None and "N" in read.cigarstring:
                spliced_read_num += 1
            total_intron += read_block.introns["intron_num"]
            if read_miss_exon_flag:
                misaligned_read_num += 1
            if read_fix_flag:
                fix_read_num += 1
        except Exception as e:
            logging.exception(e)
            error_case += 1
            error_bam.write(read)
    stat_dict = {"read_count": read_i, "spliced_read_num": spliced_read_num, "total_intron": total_intron,
                 "misaligned_read_num": misaligned_read_num,
                 "misaligned_intron_num": misaligned_intron_num,
                 "fix_read_num": fix_read_num, "fix_intron_num": fix_intron_num}
    # The read_count in stat_dict is the count of the whole bam file, not the part in the process.
    res_dict = {"filter_list": filter_list, "stat_dict": stat_dict}
    ori_bam.close()
    error_bam.close()
    if error_case == 0:
        os.remove(error_bam_path)
    if not only_region:
        realign_bam.close()
        res_dict["out_bam_path_proc"] = out_bam_path_proc
    return res_dict


def run_multiprocess(
    ori_bam_path, ref_fa_path, annot_bed_path, out_bam_path, out_realign_region_path,
    small_exon_size=80, flank_len=20, ignore_strand=True, set_tag=True, nprocess=1,
    delta_ratio_thd=0.5, simplify=True, float_flank_len=False, only_region=False, debug_mode=True
):
    start_time = time.process_time()
    start_time2 = time.time()
    annot_bed = parse_bed(annot_bed_path)
    annot_bed = annot_bed[annot_bed["exon_num"] != 1]  # Filter out singleten
    # In order to avoid repetitive computation, we extend the bed file in advance.
    annot_exon = extend_tx(annot_bed)
    read_count = int(pysam.view("-c", ori_bam_path).strip())
    filter_list = []
    if not only_region:
        out_bam_list = []
    spliced_read_num = 0
    total_intron = 0
    misaligned_read_num = 0
    misaligned_intron_num = 0
    fix_read_num = 0
    fix_intron_num = 0
    proc_list = []  # contain the list of proc_region
    per_read_num = int(read_count / nprocess)
    for i in range(nprocess):
        proc_region = {
            "proc_id": i,
            "start_read_id": i * per_read_num,
            "end_read_id": (i+1) * per_read_num
        }
        proc_list.append(proc_region)
    proc_list[-1]["end_read_id"] = read_count  # reset the end
    # def find_fix_proc(proc_region):
    #     # define a closure to use find_fix in multiprocessing.Pool instance
    #     res_dict = find_fix(
    #         proc_region, ori_bam_path, out_bam_path, ref_fa_path, annot_bed, score_matrix,
    #         small_exon_size=small_exon_size, flank_len=flank_len, ignore_strand=ignore_strand,
    #         delta_ratio_thd=delta_ratio_thd, simplify=simplify, float_flank_len=float_flank_len,
    #         only_region=only_region
    #     )
    #     return res_dict
    con_args = (
        ori_bam_path, out_bam_path, ref_fa_path, annot_bed, annot_exon,
        small_exon_size, flank_len, ignore_strand, set_tag,
        delta_ratio_thd, simplify, float_flank_len, only_region
    )
    pool = mp.Pool(processes=nprocess)
    if debug_mode:
        multi_res = [pool.apply_async(find_fix_debug, (proc_region, *con_args))
                     for proc_region in proc_list]
    else:
        multi_res = [pool.apply_async(find_fix, (proc_region, *con_args))
                     for proc_region in proc_list]
    pool.close()
    pool.join()
    res_list = [res.get() for res in multi_res]
    # res_list = pool.map(find_fix_proc, proc_list)
    for res_dict in res_list:
        filter_list.extend(res_dict["filter_list"])
        if not only_region:
            out_bam_list.append(res_dict["out_bam_path_proc"])
        spliced_read_num += res_dict["stat_dict"]["spliced_read_num"]
        total_intron += res_dict["stat_dict"]["total_intron"]
        misaligned_read_num += res_dict["stat_dict"]["misaligned_read_num"]
        misaligned_intron_num += res_dict["stat_dict"]["misaligned_intron_num"]
        fix_read_num += res_dict["stat_dict"]["fix_read_num"]
        fix_intron_num += res_dict["stat_dict"]["fix_intron_num"]
    if res_list[0]["stat_dict"]["read_count"] != read_count:
        print(
            "Warning: the total read count is not consist in pysam!",
            res_list[0]["stat_dict"]["read_count"],
            read_count
        )
    realign_df = pd.DataFrame(filter_list)
    realign_df.to_csv(out_realign_region_path, sep="\t",
                      index=False, quoting=csv.QUOTE_NONE)
    # Merge all output bam file in processes.
    pysam.merge("-cp", "-h", ori_bam_path, out_bam_path, *out_bam_list)
    for out_bam in out_bam_list:
        os.remove(out_bam)  # Remove tmp file
    end_time = time.process_time()
    end_time2 = time.time()
    print("Input BAM information:")
    print("\tTotal number of read alignment records (equal to \"samtools view -c bam\"):", read_count)
    print("\tTotal number of spliced read alignment records (containing \"N\" in CIGAR strings):", spliced_read_num)
    print("\tTotal number of intron regions on read alignment records:", total_intron)
    print("The result of finding potential misaligned regions (PMRs):")
    print("\tThe number of PMRs:", misaligned_intron_num)
    if total_intron != 0:
        print("\t\tPercentage to total number of intron regions: {0:.2f}%".format(misaligned_intron_num/total_intron*100))
    print("\tThe number of read alignment records containing PMRs:", misaligned_read_num)
    if read_count != 0:
        print("\t\tPercentage to total number of read alignment records: {0:.2f}%".format(misaligned_read_num/read_count*100))
    if spliced_read_num != 0:
        print("\t\tPercentage to total number of spliced read alignment records: {0:.2f}%".format(misaligned_read_num/spliced_read_num*100))
    print("The result of realignment (judged as misaligned cases):")
    print("\tThe number of fixed PMRs:", fix_intron_num)
    if total_intron != 0:
        print("\t\tPercentage to total number of intron regions: {0:.2f}%".format(fix_intron_num/total_intron*100))
    print("\tThe number of fixed read alignment records:", fix_read_num)
    if read_count != 0:
        print("\t\tPercentage to total number of read alignment records: {0:.2f}%".format(fix_read_num/read_count*100))
    if spliced_read_num != 0:
        print("\t\tPercentage to total number of spliced read alignment records: {0:.2f}%".format(fix_read_num/spliced_read_num*100))
    print("Main process time used: {:.2f}s".format(end_time - start_time))
    print("Time used: {:.2f}s".format(end_time2 - start_time2))
