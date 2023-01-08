"""
The script edit from finc_multi_small_exon.R in the original codes
Find the possible missed small exons, with the threshold: small_exon_size <= 100 and abs(delta_ratio) <= 0.5. 
Each realign region will be output as a dict object with the following essential key values:
    read_i(query_i), read_ith_intron(query_ith_intron), tx_i(sub_i), contain_tx_ith_exon, sum_exon_size, contain_exon_num, chrom, strand, 
    realign_start, realign_end, read_intron_start(left_query_end), read_intron_end(right_query_start), tx_leftexon_end(left_sub_end), 
    tx_rightexon_start(right_sub_start), small_exon_starts, small_exon_ends, margin_len, margin_len_mod, delta_ratio, delta_ratio_mod
Each realign region corresponds to one read intron cover one or multiple continuous small exons of one transcript isoforms.
For one read intron overlaps with multiple transcripts isoforms, we can only keep one of them for realignment using the minimum of 
    abs(delta_ratio_mod) (param: simplify = True). As delta_ratio_mod is calculated in the later step, we will use the
    simplify parameter in the fix_small_exon.py. It can reduce the time of realignment if there are many regions overlaps with multiple 
    transcript isoforms. We can also use all of the realign region for realignment (param: simplify = False). The realignment result
    will only return the alignment with best realign score.
"""


import numpy as np
import pandas as pd


EMPTY_SS = pd.Series([])


def condense_exon(exons):
    """
    The function condenses the input dataframe into a series.
    The exon_i, exon_start, exon_end will be condensed into list.
    """
    exons.reset_index(inplace=True)
    cond_exons = exons.iloc[0].drop(["exon_i", "exon_start", "exon_end"])
    cond_exons["exon_ids"] = list(exons["exon_i"])
    cond_exons["exon_starts"] = list(exons["exon_start"])
    cond_exons["exon_ends"] = list(exons["exon_end"])
    return cond_exons


def find_overlap_tx(read_block, bed_df, ignore_strand=True):
    """
    The function use read_block to find which transcripts in bed_df overlaps with the range of introns in the read.
    Return the overlapped transcript index.
    """
    bed_df = bed_df[bed_df["chrom"] == read_block.chrom]
    if not ignore_strand:
        bed_df = bed_df[bed_df["strand"] == read_block.strand]
    left_flag = bed_df["end"] > min(read_block.introns["intron_starts"])
    right_flag = bed_df["start"] < max(read_block.introns["intron_ends"])
    tx_i = bed_df[left_flag & right_flag].index
    return tx_i


def interval_overlap(i1, i2):
    """
    The function judge whether two intervals overlap with each other. i1 and i2 should be list or tuple with two elements.
    We require the first element of the interval is start and the second is end. start < end. start and end are 0-based
    Return a bool. True for overlap.
    """
    if i1[0] > i1[1] or i2[0] > i2[1]:
        raise Exception("Error: not standard interval!", i1, i2)
    # We find minimap2 might have reads with intron connect insertion and then connect intron.
    if i1[0] == i1[1] or i2[0] == i2[1]:
        flag = False
    else:
        flag = i1[1] > i2[0] and i1[0] < i2[1]
    return flag


def find_cover_exons(read_block, exons, flank_len=20, small_exon_size=100):

    global EMPTY_SS

    def merge_check(row, ith_intron, intron_start, intron_end):
        # We need to filter out the ol_exon with exons at the end and the flank exons not overlaps with read exons.
        # Check the continuity of ol_exon (as we require exon size ≤ small_exon_size).
        end_check = (not 0 in row.exon_ids) and (
            not row.exon_num-1 in row.exon_ids)
        if not end_check:
            return EMPTY_SS
        leftmost_exon_id = min(row.exon_ids)
        rightmost_exon_id = max(row.exon_ids)
        contain_exon_num = len(row.exon_ids)
        cont_check = rightmost_exon_id-leftmost_exon_id+1 == contain_exon_num
        if not cont_check:
            return EMPTY_SS
        exon_sizes = pd.Series(row.exon_ends) - pd.Series(row.exon_starts)
        size_check = all(exon_sizes <= small_exon_size)
        if not size_check:
            return EMPTY_SS
        left_tx_exon = exons.loc[(
            row.tx_i, leftmost_exon_id-1), ["exon_start", "exon_end"]]
        left_read_exon = [read_block.exons["exon_starts"]
                          [ith_intron], read_block.exons["exon_ends"][ith_intron]]
        left_check = interval_overlap(left_tx_exon, left_read_exon)
        right_tx_exon = exons.loc[(
            row.tx_i, rightmost_exon_id+1), ["exon_start", "exon_end"]]
        right_read_exon = [read_block.exons["exon_starts"]
                           [ith_intron+1], read_block.exons["exon_ends"][ith_intron+1]]
        right_check = interval_overlap(right_tx_exon, right_read_exon)
        overlap_check = left_check & right_check
        if not overlap_check:
            return EMPTY_SS
        realign_start = max(min(
            left_tx_exon[1], left_read_exon[1]) - flank_len, left_tx_exon[0], left_read_exon[0])
        realign_end = min(max(
            right_tx_exon[0], right_read_exon[0]) + flank_len, right_tx_exon[1], right_read_exon[1])
        realign_s = pd.Series({
            "read_i": read_block.read_i,
            "read_ith_intron": ith_intron,
            "tx_i": row.tx_i,
            "contain_tx_ith_exons": row.exon_ids,
            "exon_sizes": list(exon_sizes),
            "sum_exon_size": sum(exon_sizes),
            "contain_exon_num": contain_exon_num,
            "chrom": row.chrom,
            "read_strand": read_block.strand,
            "tx_strand": row.strand,
            "realign_start": realign_start,
            "realign_end": realign_end,
            "tx_leftexon_end": left_tx_exon[1],
            "tx_rightexon_start": right_tx_exon[0],
            "small_exon_starts": row.exon_starts,
            "small_exon_ends": row.exon_ends,
            "read_intron_start": intron_start,
            "read_intron_end": intron_end,
            "margin_len": intron_start - left_tx_exon[1] + right_tx_exon[0] - intron_end})
        realign_s["delta_ratio"] = (
            realign_s["margin_len"] - realign_s["sum_exon_size"]) / realign_s["sum_exon_size"]
        return realign_s
    #
    realign_list = []
    for ith_intron, intron_start, intron_end in zip(range(read_block.introns["intron_num"]),
                                                    read_block.introns["intron_starts"], read_block.introns["intron_ends"]):
        left_flag = exons["exon_start"] >= intron_start
        right_flag = exons["exon_end"] <= intron_end
        cover_exon = exons[left_flag & right_flag]
        if not cover_exon.empty:
            cond_cover_exon = cover_exon.groupby(
                level="tx_i").apply(condense_exon)
            for rt in cond_cover_exon.itertuples():
                realign_ss = merge_check(
                    rt, ith_intron, intron_start, intron_end)
                if not realign_ss.empty:
                    realign_list.append(realign_ss)
    if len(realign_list) == 0:
        return False
    else:
        return realign_list


def find_missed_exons(read_block, bed_df, exons, ignore_strand=True, flank_len=20, small_exon_size=100):
    """
    The function find transcript exons which overlap with the intron regions.
    Return a dataframe, each row is one intron with overlapped exons of one transcript.
    """
    if read_block.exons["exon_num"] == 1:  # Filter out singleten
        return False
    ol_tx_i = find_overlap_tx(read_block, bed_df, ignore_strand=ignore_strand)
    if ol_tx_i.empty:
        return False
    # Only used the exons on the transcript which overlaps with the introns to increase speed
    exons = exons.loc[(ol_tx_i, slice(None)), :]
    realign_list = find_cover_exons(
        read_block, exons, flank_len=flank_len, small_exon_size=small_exon_size)
    return realign_list
