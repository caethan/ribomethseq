"""
Take an iterator over chrom, position, and raw read-end counts and produce a second iterator 
with chrom, position, and scores, as described in Birkedal et al.
"""

from __future__ import division

import numpy as np

from ribomethseq.wiggle import WiggleReader, WiggleWriter


FLANK_WIDTH = 6
EMPTY_VALUE = 0


def split_flanks(end_counts):
    if not len(end_counts) % 2 == 1:
        raise ValueError('Full window must be odd')
    delta = len(end_counts) // 2
    left_flank = end_counts[:delta]
    center = end_counts[delta]
    right_flank = end_counts[-delta:]
    return (left_flank, 
            center, 
            right_flank)


def calculate_score_A(end_counts):
    left_flank, center, right_flank = split_flanks(end_counts)
    numerator = 2 * center + 1
    denomenator = abs(np.mean(left_flank) - np.std(left_flank, ddof=1)) / 2 + center + \
                  abs(np.mean(right_flank) - np.std(right_flank, ddof=1)) / 2 + 1
    return max(1 - numerator / denomenator, 0.0)


def count_wig_to_score_wig(count_file, score_file, score_func):
    writer = None
    for chrom, position, score in iterate_scores_from_wiggle(count_file, score_func):
        if writer is None:
            writer = WiggleWriter(chrom, 1, score_file)
            writer.write_header()
        writer.write_score(position + 1, score) # Switch back to 1-based wiggle filename
    writer.outfile.close()

    
    
def iterate_scores_from_wiggle(count_file, score_func):
    for chrom, position, end_counts in iterate_windows_from_wiggle(count_file):
        yield chrom, position, score_func(end_counts)

    
    
def iterate_windows_from_wiggle(filename):
    reader = WiggleReader(filename)
    data = [(int(position), float(score)) for _, position, score in reader.iterate_wiggle(True)]
    
    def naive_window(index):
        center_position, score = data[index]
        output = [score]
        # Read back
        left_flank = []
        test_index = index - 1
        for i in range(FLANK_WIDTH):
            expected_position = center_position - (i + 1)
            while True:
                if test_index < 0:
                    new_score = EMPTY_VALUE
                    break
                else:
                    test_pos, test_score = data[test_index]
                    if test_pos == expected_position:
                        new_score = test_score
                        break
                    elif test_pos > expected_position:
                        test_index -= 1
                        continue
                    else:
                        new_score = EMPTY_VALUE
                        break
            left_flank.append(new_score)
        left_flank = left_flank[::-1]
        # Read forward
        right_flank = []
        test_index = index + 1
        for i in range(FLANK_WIDTH):
            expected_position = center_position + i + 1
            while True:
                if test_index >= len(data):
                    new_score = EMPTY_VALUE
                    break
                else:
                    test_pos, test_score = data[test_index]
                    if test_pos == expected_position:
                        new_score = test_score
                        break
                    elif test_pos < expected_position:
                        test_index += 1
                        continue
                    else:
                        new_score = EMPTY_VALUE
                        break
            right_flank.append(new_score)
        return left_flank + output + right_flank
        
    for index in range(len(data)):
        yield reader.current_chrom, data[index][0], naive_window(index)
            
        