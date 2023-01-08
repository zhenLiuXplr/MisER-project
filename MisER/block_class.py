import re
import warnings

import numpy as np

from copy import deepcopy


"""
The Cigar class is try to parse cigar string and achieve fast indexing.
    We define three kinds of positions of an aligned read:
    1. reference position: The position based on the reference. Normally, start or end without any suffix is reference position.
    2. read position (suffix: rdp or rdps): The position based on the read, which means the ith of the read base at this position.
    3. extend cigar position (suffix: ecp): The position based on the extend cigar.
    4. cigar index (suffix: ci): The index of cigar list.
    All of these positions are 0-based.
    Examples:
            * 
        012345678 (ecp)
        SSMMDDMMM
            *
          0123456 (pos)
        --ATTTGCC
              *
        0123  456 (rdp)
        GGAT--GGC

        If the alignment come the first star site (ecp=4, pos=2, rdp=4), The cigar character 'D' at ecp=4 will determine the
        alignment stage of 'T' at pos=2 and 'G' at ecp=4. As 'D' means deletion, so ecp and pos increase one and 
        point to the next one, while rdp remain unchanged rdp=4 .

               *
        012345678 (ecp)
        SSMMDDMMM
               *
          0123456 (pos)
        --ATTTGCC
               *
        0123  456 (rdp)
        GGAT--GGC

        If the alignment come the second star site (ecp=7, pos=5, rdp=5), The cigar character 'M' at ecp=7 will determine the
        alignment stage of 'C' at pos=5 and 'G' at ecp=5. As 'M' means match, so ecp, pos and rdp all increase one and 
        point to the next one.

    The relationship between pos, rdp, ecp can be writed as functions:
    pos = f1(ci)
    rdp = f2(ci)
    ci = g(pos, rdp)
"""


class CigarSliceError(Exception):
    """
    The class is used to raise an error when the cigar slice operation overstepped or start_ci > end_ci
    """

    def __init__(self, start_ci, end_ci, start_bound, end_bound):
        self.start_ci = start_ci
        self.end_ci = end_ci
        self.start_bound = start_bound
        self.end_bound = end_bound

    def __str__(self):
        return "cigar slice operation overstepped or start_ci > end_ci. \
                start_ci:{0}, start_bound:{1}, end_ci:{2}, end_bound:{3}" \
                    .format(self.start_ci, self.start_bound, self.end_ci, self.end_bound)

class CigarDefectError(Exception):
    def __init__(self, qname):
        self.qname = qname
    
    def __str__(self):
        return "Segement:{0} missing cigar field in the BAM file.".format(self.qname)


class Cigar(object):
    """
    We convert cigar string to cigar list. e.g. '6S55M77N' -> [[6, 'S'], [55, 'M'], [77, 'N']]
    In order to accelerate the tracing process. We abandon the method to use the extend cigar string and 
    the absolute extend cigar position (ecp). Instead, we index the cigar list at two levels (i, j).
    The first index i means the index of cigar list. e.g. i=1 trace to character group [55, 'M'] in the cigar list above
    The second index j means the index of  character. e.g. i=1 and j=10 trace to the 11th 'M' in [55, 'M']
    The cigar index (ci) is a tuple (i, j).
    """
    regexNum = re.compile('[0-9]+')
    regexChr = re.compile('[A-Z=]')
    cigarSet = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'}

    def __init__(self, cstr):
        num_list = [int(x) for x in Cigar.regexNum.findall(cstr)]
        chr_list = Cigar.regexChr.findall(cstr)
        if not(set(chr_list) <= Cigar.cigarSet):
            raise Exception(
                "error: CIGAR contain the pattern exceed the assumption!", set(chr_list))
        self.str = cstr
        self.list = [[x, y] for x, y in zip(num_list, chr_list)]

    def __add_cumsum(self):
        self.cumsum_list = np.cumsum([x[0] for x in self.list])
    
    def getstr(self):
        self.str = ''
        for n, c in self.list:
            self.str = self.str + str(n) + c
        return self.str
    
    def ecp2ci(self, ecp):
        """
        Convert absolute extend cigar position to cigar index.
        Slow.
        """
        if not hasattr(self, 'cumsum_list'):
            self.__add_cumsum()
        if ecp < 0 or ecp > self.cumsum_list[-1]:
            raise Exception("Error: ecp out of bond!", ecp, self.cumsum_list[-1])
        for i, csum in zip(range(len(self.cumsum_list)), self.cumsum_list):
            if ecp < csum:
                break
        j = ecp - (self.cumsum_list[i-1] if i > 0 else 0)
        return (i, j)

    def ci2ecp(self, ci):
        """
        Convert cigar index to absolute extend cigar position.
        Quick.
        """
        if not hasattr(self, 'cumsum_list'):
            self.__add_cumsum()
        ecp = self.cumsum_list[ci[0] - 1] if ci[0] > 0 else 0
        ecp += ci[1]
        return ecp
    
    def standardize_list(cigar_list):
        """
        Check and standardize the cigar list. Include the following operation:
        1. merge adjacent cigar group with same kind of character.
        2. check and ensure the cigar number is positive. Delete 0 and raise error if negative.
        """
        if len(cigar_list) == 0:
            return []
        if cigar_list[0][0] == 0:
            return []
        standard_list = [deepcopy(cigar_list[0])]
        for n, c in cigar_list[1:]:
            if n == 0:
                continue
            elif n < 0:
                raise Exception("Error: negative cigar number!", cigar_list)
            if c == standard_list[-1][1]:
                standard_list[-1][0] += n
            else:
                standard_list.append([n, c])
        return standard_list

    def standardize(self):
        self.list = Cigar.standardize_list(self.list)
                

    def cslice(self, start_ci=None, end_ci=None):
        """
        start_ci: start cigar index
        end_ci: end cigar index
        Return a slice of cigar the list.
        Note: the cigar index is 0-based and the slice operation is half-open (i is closed but j is left-closed and right-open).
        An example:
        tc = Cigar("6S7M3D8M")
        tc.list:
        [[6, 'S'], [7, 'M'], [3, 'D'], [8, 'M']]
        tc.cslice(end_ci = (2, 0)) will return [[6, 'S'], [7, 'M']] 
        ATTENTION: The index (1, 7) is deprecated, we only allow the last cigar group can have index over-bound.
        """
        list_len = len(self.list)
        if start_ci is None:
            start_ci = (0, 0)
        if end_ci is None:
            end_ci = (list_len-1, self.list[-1][0])
        if start_ci == (list_len-1, self.list[-1][0]):
            return []
        # ci check
        if end_ci[0] == list_len-1:
            if not(0 <= start_ci[0] < list_len) or not(0 <= start_ci[1] < self.list[start_ci[0]][0]) \
            or not(0 <= end_ci[1] <= self.list[end_ci[0]][0]) or start_ci[0] > end_ci[0]:
                raise CigarSliceError(start_ci, end_ci, self.list[start_ci[0]][0], self.list[end_ci[0]][0])
        else:
            if not(0 <= start_ci[0] < list_len) or not(0 <= start_ci[1] < self.list[start_ci[0]][0]) \
            or not(end_ci[0] < list_len) or not(0 <= end_ci[1] < self.list[end_ci[0]][0]) or start_ci[0] > end_ci[0]:
                raise CigarSliceError(start_ci, end_ci, self.list[start_ci[0]][0], self.list[end_ci[0]][0])
        slice_list = deepcopy(self.list[start_ci[0]:end_ci[0]+1])
        if start_ci[0] == end_ci[0]:  # slice on the same one character group
            if end_ci[1] < start_ci[1]:
                raise CigarSliceError(start_ci, end_ci, self.list[start_ci[0]][0], self.list[end_ci[0]][0])
            slice_list[0][0] = end_ci[1] - start_ci[1]
        else:
            slice_list[0][0] = slice_list[0][0] - start_ci[1]
            if end_ci[1] == 0:
                slice_list = slice_list[:-1]
            else:
                slice_list[-1][0] = end_ci[1]
        slice_list = Cigar.standardize_list(slice_list)
        return slice_list

    def trace_ci_rdp(self, start_ci, start_pos, end_pos, start_rdp, end_rdp=None):
        """
        The defination of pos, rdp see The class ReadBlock.
        Actuallly, after the initial value start_ci, start_pos, start_rdp setted, we still need two arguments to 
        calculate end_ci (ci = g(pos, rdp)), end_pos and end_rdp. But usually, we only provide one argument, as 
        we only use end_pos to calculate end_ci, we will return end_ci and end_rdp as the pos first reach end_pos. 
        Return (end_ci, end_rdp, reach_flag)
        """
        pos = start_pos
        rdp = start_rdp
        ecp = self.ci2ecp(start_ci)
        reach_flag = False
        if start_pos == end_pos:
            reach_flag = True
            return (start_ci, start_rdp, reach_flag)
        if start_pos < end_pos:
            # forward tracing
            cigar_slice = self.cslice(start_ci=start_ci)
            for n, c in cigar_slice:
                if pos == end_pos:
                    reach_flag = True
                    break
                if c in {'M', 'X', '='}:
                    for j in range(n):
                        pos += 1
                        rdp += 1
                        ecp += 1
                        if pos == end_pos:
                            reach_flag = True
                            break
                elif c in {'S', 'I'}:
                    rdp += n
                    ecp += n
                elif c in {'D', 'N'}:
                    for j in range(n):
                        pos += 1
                        ecp += 1
                        if pos == end_pos:
                            reach_flag = True
                            break
                else:
                    ecp += n
        else:
            # backward tracing
            cigar_slice = self.cslice(end_ci=start_ci)
            for n, c in reversed(cigar_slice):
                if pos == end_pos:
                    reach_flag = True
                    break
                if c in {'M', 'X', '='}:
                    for j in range(n):
                        pos -= 1
                        rdp -= 1
                        ecp -= 1
                        if pos == end_pos:
                            reach_flag = True
                            break
                elif c in {'S', 'I'}:
                    for j in range(n):
                        rdp -= n
                        ecp -= n
                elif c in {'D', 'N'}:
                    for j in range(n):
                        pos -= 1
                        ecp -= 1
                        if pos == end_pos:
                            reach_flag = True
                            break
                else:
                    ecp -= n
        end_ci = self.ecp2ci(ecp)
        end_rdp = rdp
        return (end_ci, end_rdp, reach_flag)


class Cigar_fromlist(Cigar):
    """
    Create cigar object from cigar list rather than string.
    """

    def __init__(self, cigar_list):
        self.list = cigar_list

    

class ReadBlock(object):
    """
    The class is try to parse one read from BAM file (pysam.AlignedSegment) and return attribute of the blocks.
    To accelerate tracing process, we also return cigar list index to trace the introns.
    """

    def __init__(self, read, read_i):
        self.read_i = read_i
        self.chrom = read.reference_name
        self.strand = '-' if read.is_reverse else '+'
        self.cigar = Cigar(read.cigarstring)
        pos = read.reference_start
        rdp = 0
        self.starts = [pos]
        self.ends = []
        for n, c in self.cigar.list:
            if c in {'M', '=', 'X'}:
                pos += n
                rdp += n
            elif c in {'S', 'I'}:
                rdp += n
            elif c in {'D', 'N'}:
                if c == 'N':
                    self.ends.append(pos)
                    pos += n
                    self.starts.append(pos)
                else:
                    pos += n
        self.ends.append(pos)
        self.exons = {"exon_num": len(
            self.starts), "exon_starts": self.starts, "exon_ends": self.ends}
        self.introns = {"intron_num": len(self.starts)-1, "intron_starts": self.ends[:-1], "intron_ends": self.starts[1:]}
