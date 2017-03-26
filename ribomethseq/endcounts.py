"""Tools for counting up the end positions of reads in an input BAM file.

Some tricky bits:
    - Input BAM files and output WIG files use different coordinate systems (0 based for BAM,
      1-based for WIG).
    - The 5' end of reads needs to be decremented by 1 due to the chemistry of the
      methylation/fragmentation
    - We assume paired reads, so we only read out the 5' end of upstream reads and
      the 3' end of downstream reads.  At some point we may wish to add an option to 
      correctly handle unpaired reads.  
"""

import argparse
import pysam
import numpy as np

from ribomethseq import wiggle


def is_good_paired_read(read):
    return all([read.is_proper_pair and
                not read.is_qcfail and
                not read.is_unmapped and
                not read.mate_is_unmapped and
                not read.is_supplementary and
                not read.is_duplicate])

def is_proper_pair(read):
    return (read.is_read1 and not read.is_reverse) or (read.is_read2 and read.is_reverse):


def get_read_coords(read):
    chrom = read.reference_name
    # We're converting from the 0-based half open BED file to the 1-based closed WIG file 
    # format here, incrementing the start coordinate by 1
    start = read.reference_start + 1
    stop = read.reference_end
    return chrom, start, stop


def create_iterator(chrom, end_counter):
    for position in range(max(end_counter.keys()) + 1):
        yield chrom, position, end_counter[position]


def get_endcount_iterators(in_bam_file):
    # make sure the input file is a BAM and not a SAM
    if not in_bam_file.endswith('.bam'):
        raise ValueError('Must provide a BAM file input')
    ends = {}
    samfile = pysam.AlignmentFile(in_bam_file, "rb")
    # Without args, this gets all mapped reads
    for read in samfile.fetch():
        if not is_good_paired_read or not is_proper_pair(read):
            continue
    
        read_type = 'read1' if read.is_read1 else 'read2'
        chrom, start, stop = get_read_coords(read)
        # Downshift the 5' position by 1 because of the fragmentation chemistry
        start -= 1
    
        if not chrom in ends:
            ends[chrom] = {'5-prime': Counter(),
                           '3-prime': Counter(),
                           'total': Counter(),}
        if read_type == 'read1':
            ends[chrom]['5-prime'][start] += 1
            ends[chrom]['total'][start] += 1
        elif read_type == 'read2':
            ends[chrom]['3-prime'][stop] += 1
            ends[chrom]['total'][stop] += 1
            
    iterators = {}
    for chrom in ends.keys():
        iterators[chrom] = (
            create_iterator(chrom, ends[chrom]['5-prime']),
            create_iterator(chrom, ends[chrom]['3-prime']),
            create_iterator(chrom, ends[chrom]['total'])
        )
    return iterators