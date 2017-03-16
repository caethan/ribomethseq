"""
Tools for reading and writing the Wiggle Track Format:
https://genome.ucsc.edu/goldenpath/help/wiggle.html

Wiggle format positions are 1-relative.

Returns positions that match the BED format, which is zero-based, half-open.

Currently only supports variableStep format.
"""


MODE_VARIABLE = 'variable'


class WiggleReader:
    
    def __init__(self, filename):
        self.filename = filename
        self.current_chrom = None
        self.current_pos = -1
        self.current_span = -1
        self.current_step = -1
        self.mode = None
        
    
    def parse_header(self, line):
        return dict(field.split('=') for field in line.split()[1:])    

    
    def iterate_wiggle(self, iterate_by_position=True):
        with open(self.filename, 'r') as infile:
            for line in infile:
                if line.isspace() or line.startswith('#'):
                    continue
                if line.startswith('track') or line.startswith('browser'):
                    continue
                # Read and set mode
                if line.startswith('variableStep'):
                    if not self.mode is None:
                        raise ValueError('Cannot set mode more than once per file')
                    header = self.parse_header(line)
                    self.mode = MODE_VARIABLE
                    self.current_chrom = header['chrom']
                    self.current_span = int(header.get('span', 1))
                    continue
                # Read data with the active mode
                if self.mode == MODE_VARIABLE:
                    start, score = line.split()
                    start = int(start) - 1 # Again, switch from 1-based to 0-based
                    score = float(score)
                    if iterate_by_position:
                        for position in range(start, start + self.current_span):
                            yield self.current_chrom, position, score
                    else:
                        yield self.current_chrom, start, start + self.current_span, score
                else:
                    raise ValueError('Unexpected line: {}'.format(line))
                    

class WiggleWriter:
    def __init__(self, chrom, span, filename):
        self.span = span
        self.chrom = chrom
        self.mode = MODE_VARIABLE
        self.outfile = open(filename, 'w')
        
    
    def write_header(self):
        self.outfile.write('variableStep chrom={chrom} span={span}\n'.format(chrom=self.chrom, span=self.span))

    
    def write_score(self, position, score, round=False):
        if round:
            self.outfile.write('{position}\t{score:.4f}\n'.format(**locals()))
        else:
            self.outfile.write('{position}\t{score}\n'.format(**locals()))
        

    def write_from_iterator(self, filename, iterator):
        self.write_header()
        for locus in iterator:
            if len(locus == 4):
                chrom, start, stop, score = locus
                span = stop - start
            else:
                chrom, start, score = locus
                span = 1
            if not span == self.span:
                raise ValueError('Span mismatch: {}, {}'.format(span, self.span))
            if not chrom == self.chrom:
                raise ValueError('Chrom mismatch: {}, {}'.format(chrom, self.chrom))
            self.write_score(start, score)