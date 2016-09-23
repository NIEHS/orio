import os

from .base import Validator


class AnalysisValidator(Validator):

    ANCHOR_OPTIONS = ('start', 'end', 'center')
    MIN_BIN_NUMBER = 1
    MAX_BIN_NUMBER = 250
    MIN_BIN_SIZE = 1
    MAX_WINDOW_SIZE = 100000

    # Possible 'empty' names used in a BED file
    DUMMY_NAMES = ('.', )

    def __init__(self, bin_anchor, bin_start, bin_number,
                 bin_size, feature_bed, chrom_sizes,
                 stranded_bed):
        super().__init__()

        # Type checks
        assert isinstance(bin_anchor, str)
        assert isinstance(bin_start, int)
        assert isinstance(bin_number, int)
        assert isinstance(bin_size, int)
        assert os.path.exists(feature_bed)
        assert os.path.exists(chrom_sizes)

        # Set instance variables
        self.bin_anchor = bin_anchor
        self.bin_start = bin_start
        self.bin_number = bin_number
        self.bin_size = bin_size
        self.feature_bed_fn = feature_bed
        self.chrom_sizes_fn = chrom_sizes
        self.stranded_bed = stranded_bed

    def validate(self):
        self.check_input_domains()
        self.check_if_outside()

    def check_input_domains(self):
        if self.bin_size < self.MIN_BIN_SIZE:
            self.add_error('Bin size {} less than minimum {}'.format(
                self.bin_size, self.MIN_BIN_SIZE))

        window_size = self.bin_size * self.bin_number
        if window_size > self.MAX_WINDOW_SIZE:
            self.add_error('Window size {} exceeds maximum {}'.format(
                window_size, self.MAX_WINDOW_SIZE))

        if self.bin_number < self.MIN_BIN_NUMBER:
            self.add_error('Bin number {} less than minimum {}'.format(
                self.bin_number, self.MIN_BIN_NUMBER))

        if self.bin_number > self.MAX_BIN_NUMBER:
            self.add_error('Bin number {} exceeds maximum {}'.format(
                self.bin_number, self.MAX_BIN_NUMBER))

    def read_chrom(self, chrom_sizes_fn):
        chrom_sizes = dict()
        with open(chrom_sizes_fn) as f:
            for line in f:
                chromosome, size = line.strip().split('\t')
                chrom_sizes[chromosome] = int(size)
        return chrom_sizes

    def check_header(self, line):
        # Check to see if line is header
        if line == '\n':
            return True
        elif line[0] == '#':
            return True
        elif line.split()[0].lower() in ('track', 'browser'):
            return True
        else:
            return False

    def check_if_outside(self):
        outside_message = 'Feature window extends off chromosome. BED entry will be removed from analysis.'  # noqa

        window_size = self.bin_size * self.bin_number
        chrom_sizes = self.read_chrom(self.chrom_sizes_fn)
        entry_num = 0
        with open(self.feature_bed_fn) as f:
            for line in f:
                if not self.check_header(line):
                    bed_fields = len(line.strip().split())
                    chromosome, start, end = line.strip().split()[0:3]
                    start = int(start) + 1  # Convert from 0-based to 1-based
                    end = int(end)
                    center = int((start + end) / 2)

                    entry_name = None
                    if bed_fields >= 4:
                        name = line.strip().split()[3]
                        if name not in self.DUMMY_NAMES:
                            entry_name = name

                    if self.stranded_bed:
                        if bed_fields >= 6:  # Contains strand information?
                            strand = line.strip().split()[5]
                            if not (strand == '+' or strand == '-'):
                                self.add_error('Strand column is not valid.')
                        else:
                            self.add_error('BED file lacks strand column')
                    else:
                        strand = 'AMBIG'

                    # Define start and end points for window
                    if self.bin_anchor == 'center':
                        if strand == '+' or strand == 'AMBIG':
                            window_start = center + self.bin_start
                        if strand == '-':
                            window_start = center - self.bin_start
                    elif self.bin_anchor == 'start':
                        if strand == '+' or strand == 'AMBIG':
                            window_start = start + self.bin_start
                        if strand == '-':
                            window_start = end - self.bin_start
                    elif self.bin_anchor == 'end':
                        if strand == '+' or strand == 'AMBIG':
                            window_start = end + self.bin_start
                        if strand == '-':
                            window_start = start - self.bin_start

                    if strand == '+' or strand == 'AMBIG':
                        window_end = window_start + window_size
                    elif strand == '-':
                        window_end = window_start - window_size

                    chrome_size = chrom_sizes.get(chromosome)
                    if chrome_size is None:
                        self.add_error('Chromosome not found {}'.format(chromosome))
                    elif strand == '+' or strand == 'AMBIG':
                        if window_start < 1 or window_end > chrome_size:
                            if entry_name:
                                self.add_warning('{}: {}'.format(
                                    entry_name, outside_message))
                            else:
                                self.add_warning('Entry number {}: {}'.format(
                                    str(entry_num), outside_message))
                    elif strand == '-':
                        if window_end < 1 or window_start > chrome_size:
                            if entry_name:
                                self.add_warning('{}: {}'.format(
                                    entry_name, outside_message))
                            else:
                                self.add_warning('Entry number {}: {}'.format(
                                    str(entry_num), outside_message))
                    entry_num += 1
