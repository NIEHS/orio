#!/usr/bin/env python

import click
import os
import tempfile
from subprocess import call
from six import string_types
import warnings

from . import utils


class BedMatrix(object):

    ANCHOR_OPTIONS = ('start', 'end', 'center')
    bigWigAverageOverBed_path = os.path.join(
        utils.get_bin_path(), "bigWigAverageOverBed")
    DUMMY_VALUES = ['.']

    def __init__(self, bigwigs, feature_bed, output_matrix, anchor, bin_start,
                 bin_number, bin_size, opposite_strand_fn, stranded_bigwigs,
                 stranded_bed, chrom_sizes, temp_dir=None):

        # Set instance variables
        self.feature_bed = feature_bed
        self.output_matrix = output_matrix
        self.anchor = anchor
        self.bin_start = bin_start
        self.bin_number = bin_number
        self.bin_size = bin_size
        self.opposite_strand_fn = opposite_strand_fn
        self.stranded_bigwigs = stranded_bigwigs
        if self.stranded_bigwigs:
            if len(bigwigs) != 2:
                raise ValueError("Two paired bigwig files are expected!!")
            self.plus_bigwig = bigwigs[0]
            self.minus_bigwig = bigwigs[1]
        else:
            if len(bigwigs) != 1:
                raise ValueError("One unstranded bigwig file is expected!!")
            self.unstranded_bigwig = bigwigs[0]
        self.stranded_bed = stranded_bed
        self.chrom_sizes = chrom_sizes
        self.temp_dir = temp_dir

        # Type checks
        assert anchor in self.ANCHOR_OPTIONS
        assert isinstance(self.bin_start, int)
        assert isinstance(self.bin_number, int)
        assert isinstance(self.bin_size, int)
        if self.opposite_strand_fn is not None:
            assert isinstance(self.opposite_strand_fn, string_types)
            if os.path.dirname(self.output_matrix) != "":
                assert os.path.exists(os.path.dirname(self.opposite_strand_fn))
            if not self.stranded_bed:
                raise ValueError("Cannot report opposite strand coverage \
                                 without stranded_bed flag!!")
        if self.stranded_bigwigs:
            assert os.path.exists(self.plus_bigwig)
            assert os.path.exists(self.minus_bigwig)
        else:
            assert os.path.exists(self.unstranded_bigwig)
        if self.chrom_sizes:
            assert os.path.exists(self.chrom_sizes)
        assert os.path.exists(self.feature_bed)
        if os.path.dirname(self.output_matrix) != "":
            assert os.path.exists(os.path.dirname(self.output_matrix))

        if self.temp_dir is None:
            self.temp_dir = tempfile.gettempdir()
        os.makedirs(self.temp_dir, exist_ok=True)

        self.execute()

    def execute(self):
        try:
            self.create_temps()
            self.make_bed()
            self.bigwig_average_over_bed()
            self.make_matrix()
        finally:
            self.cleanup()

    def tryDelete(self, fn):
        try:
            os.remove(fn)
        except FileNotFoundError:
            pass

    def create_temps(self):
        path = self.temp_dir
        self._plus_filename = tempfile.mkstemp(dir=path)[1]
        self._minus_filename = tempfile.mkstemp(dir=path)[1]
        self._unstranded_filename = tempfile.mkstemp(dir=path)[1]
        self._bedfile = tempfile.NamedTemporaryFile(mode='w', dir=path)

    def cleanup(self):
        self.tryDelete(self._plus_filename)
        self.tryDelete(self._minus_filename)
        self.tryDelete(self._unstranded_filename)
        self._bedfile.close()

    def checkInt(self, num):
        # Return int if possible, else return float
        try:
            return int(num)
        except ValueError:
            return float(num)

    @staticmethod
    def countValidBedLines(input_file):
        # Count valid feature lines in bed file
        line_count = 0
        with open(input_file) as f:
            for line in f:
                if line == "\n":
                    pass
                elif line[0] == "#":
                    pass
                elif line.split()[0].lower() in ("track", "browser"):
                    pass
                else:
                    line_count += 1
        return line_count

    @staticmethod
    def generateFeatureName(feature_header, feature_count, num_lines):
        # Generate row/feature names for bed entries
        feature_name = feature_header + "_"
        for i in range(len(str(feature_count)), len(str(num_lines))):
            feature_name = feature_name + "0"
        feature_name = feature_name + str(feature_count)
        return feature_name

    @staticmethod
    def checkHeader(line):
        # Check to see if line is header
        if line == "\n":
            return True
        elif line[0] == "#":
            return True
        elif line.split()[0].lower() in ("track", "browser"):
            return True
        else:
            return False

    @staticmethod
    def readTabName(tab_name):
        # Read feature name from tab-identifier
        feature_name = tab_name.split("_")
        del feature_name[-1]
        feature_name = "_".join(feature_name)
        return feature_name

    @staticmethod
    def readFeatureIndex(tab_name):
        return int(tab_name.split('_')[-1])

    def readChromSizes(self, in_file):
        sizes = dict()

        with open(in_file) as f:
            for line in f:
                chrom, size = line.strip().split()
                sizes[chrom] = int(size)

        return sizes

    def make_bed(self):
        """
        Open bed file, create a bed of bins.
        Create a dictionary with feature information.
        Create list with order of features in bed file.
        """
        sizes = None
        if self.chrom_sizes:
            sizes = self.readChromSizes(self.chrom_sizes)

        input_file = self.feature_bed
        feature_dict = {}
        feature_list = []
        with open(input_file, 'r') as f, \
                self._bedfile.file as OUTPUT:
            total_valid_lines = self.countValidBedLines(input_file)
            count = 0
            for line in f:
                if not self.checkHeader(line):
                    skip = False

                    # Read fields from bed
                    bed_fields = len(line.strip().split())
                    chromosome, start, end = line.strip().split()[0:3]
                    start = int(start) + 1  # 0-based to 1-based
                    end = int(end)
                    center = int((start + end) / 2)

                    name = None
                    if bed_fields >= 4:  # Contains name information?
                        name = line.strip().split()[3]
                    if name is None or name in self.DUMMY_VALUES:
                        name = self.generateFeatureName(
                            "feature", count, total_valid_lines)
                    count += 1

                    if self.stranded_bed:
                        if bed_fields >= 6:  # Contains strand information?
                            strand = line.strip().split()[5]
                            if strand in self.DUMMY_VALUES:
                                raise ValueError('Strand column is not \
                                    informative!!')
                        else:
                            raise ValueError("BED file lacks strand column!!")
                    else:
                        strand = "AMBIG"

                    # Define anchor point for window
                    if self.anchor == "center":
                        start = center
                    elif self.anchor == "start":
                        if strand == "-":
                            start = end
                    elif self.anchor == "end":
                        if strand == "+" or strand == "AMBIG":
                            start = end

                    # Create bed with bins for the given line
                    if sizes:
                        if chromosome in sizes:
                            if strand == "+" or strand == "AMBIG":
                                s = start + self.bin_start
                                e = start + self.bin_number * self.bin_size \
                                    - 1

                                if s < 1 or e > sizes[chromosome]:
                                    skip = True

                                    message = 'Window outside chromosome:\n' +\
                                        line.strip()
                                    warnings.warn(message)

                            elif strand == "-":
                                s = start - self.bin_start
                                e = start - self.bin_number * self.bin_size \
                                    + 1

                                if e < 1 or s >= sizes[chromosome]:
                                    skip = True

                                    message = 'Window outside chromosome:\n' +\
                                        line.strip()
                                    warnings.warn(message)

                        else:
                            skip = True

                            message = 'Chromosome not in sizes file:\n' + \
                                line.strip()
                            warnings.warn(message)

                    if not skip:
                        if strand == "+" or strand == "AMBIG":
                            start = start + self.bin_start
                            for i in range(self.bin_number):
                                str_ = '{}\t{}\t{}\t{}_{}\n'.format(
                                    chromosome,
                                    start - 1,
                                    start + self.bin_size - 1,
                                    name,
                                    i
                                )
                                OUTPUT.write(str_)
                                start += self.bin_size

                        elif strand == "-":
                            start = start - self.bin_start
                            for i in range(self.bin_number):
                                str_ = '{}\t{}\t{}\t{}_{}\n'.format(
                                    chromosome,
                                    start - self.bin_size,
                                    start,
                                    name,
                                    i
                                )
                                OUTPUT.write(str_)
                                start -= self.bin_size

                        feature_list.append(name)

                        feature_dict[name] = {
                            "chromosome": chromosome,
                            "start": start,
                            "end": end,
                            "strand": strand
                        }

        self.feature_order = feature_list
        self.feature_info = feature_dict

    def bigwig_average_over_bed(self):
        # Run bigWigAverageOverBed considering bed of bins
        if self.stranded_bigwigs:
            call([
                self.bigWigAverageOverBed_path,
                self.plus_bigwig,
                self._bedfile.name,
                self._plus_filename,
            ])
            call([
                self.bigWigAverageOverBed_path,
                self.minus_bigwig,
                self._bedfile.name,
                self._minus_filename,
            ])
        else:
            call([
                self.bigWigAverageOverBed_path,
                self.unstranded_bigwig,
                self._bedfile.name,
                self._unstranded_filename,
            ])

    def make_matrix(self):
        # Create output matrix files
        if self.stranded_bigwigs:
            self.make_stranded_matrix()
        else:
            self.make_unstranded_matrix()

    def readIntersectionFile(self, input_file):
        row_dict = dict()

        for feature in self.feature_order:
            row_dict[feature] = []
            for i in range(self.bin_number):
                row_dict[feature].append('0')

        for line in input_file:
            tab_name, size, covered, bed_sum, bed_mean_zero, bed_mean = \
                line.strip().split()
            bed_sum = self.checkInt(bed_sum)
            feature_name = self.readTabName(tab_name)
            feature_index = self.readFeatureIndex(tab_name)
            row_dict[feature_name][feature_index] = str(abs(bed_sum))

        return row_dict

    def parseIntersectionAsLines(self, intersection_dict):
        lines = dict()
        for feature in intersection_dict:
            lines[feature] = ''
            for i in range(len(intersection_dict[feature])):
                lines[feature] += '\t' + intersection_dict[feature][i]
        return lines

    def make_unstranded_matrix(self):
        # Create output matrix for unstranded bigwig
        with open(self._unstranded_filename, 'r') as f, \
                open(self.output_matrix, 'w') as OUTPUT:

            # Make header, write to output
            for i in range(self.bin_number):
                str_ = "\t{}:{}".format(
                    self.bin_start+i*self.bin_size,
                    self.bin_start+(i+1)*self.bin_size-1
                )
                OUTPUT.write(str_)
            OUTPUT.write("\n")

            # Make features, add to dictionary
            row_dict = self.readIntersectionFile(f)
            lines = self.parseIntersectionAsLines(row_dict)

            # Write to output based on order in feature list
            for feature_name in self.feature_order:
                OUTPUT.write(feature_name + lines[feature_name] + "\n")

    def make_stranded_matrix(self):
        """
        Create output matrix for stranded bigwig
        """
        # If opposite_strand filename exists, open opposite strand for output
        OPPOSITE = None
        if self.opposite_strand_fn:
            OPPOSITE = open(self.opposite_strand_fn, "w")

        with open(self._plus_filename) as plus_file, \
                open(self._minus_filename) as minus_file, \
                open(self.output_matrix, "w") as OUTPUT:

            # Make header
            for i in range(self.bin_number):
                OUTPUT.write("\t{}:{}".format(
                    self.bin_start+i*self.bin_size,
                    self.bin_start+(i+1)*self.bin_size-1
                ))
                if OPPOSITE:
                    OPPOSITE.write("\t{}:{}".format(
                        self.bin_start+i*self.bin_size,
                        self.bin_start+(i+1)*self.bin_size-1
                    ))

            OUTPUT.write("\n")
            if OPPOSITE:
                OPPOSITE.write("\n")

            plus_row_dict = self.readIntersectionFile(plus_file)
            minus_row_dict = self.readIntersectionFile(minus_file)

            plus_lines = self.parseIntersectionAsLines(plus_row_dict)
            minus_lines = self.parseIntersectionAsLines(minus_row_dict)

            if not self.stranded_bed:
                combined_lines = dict()
                for feature in plus_row_dict:
                    combined_lines[feature] = ''
                    for i in range(len(plus_row_dict[feature])):
                        combined_lines[feature] += '\t' + str(self.checkInt(
                                float(plus_row_dict[feature][i]) +
                                float(minus_row_dict[feature][i])))

            for feature in self.feature_order:
                OUTPUT.write(feature)
                if OPPOSITE:
                    OPPOSITE.write(feature)
                if self.stranded_bed:
                    if self.feature_info[feature]["strand"] == "+":
                        OUTPUT.write(plus_lines[feature] + '\n')
                        if OPPOSITE:
                            OPPOSITE.write(minus_lines[feature] + '\n')
                    elif self.feature_info[feature]["strand"] == "-":
                        OUTPUT.write(minus_lines[feature] + '\n')
                        if OPPOSITE:
                            OPPOSITE.write(plus_lines[feature] + '\n')
                else:
                    OUTPUT.write(combined_lines[feature] + '\n')

        if OPPOSITE:
            OPPOSITE.close()


@click.command()
@click.argument('bigwigs', nargs=-1)
@click.argument('feature_bed')
@click.argument('output_matrix')
@click.option('--temp_dir', default=None, help='Directory for temporary files',
              type=str)
@click.option('-a', '--anchor', default='center', help='Bin anchor',
              type=click.Choice(BedMatrix.ANCHOR_OPTIONS))
@click.option('-b', '--bin_start', default=-2500, help='Relative bin start',
              type=int)
@click.option('-n', '--bin_number', default=50, help='Number of bins',
              type=int)
@click.option('-s', '--bin_size', default=100, help='Size of bins',
              type=int)
@click.option('--opposite_strand_fn', type=str,
              help='Output filename for opposite strand coverage')
@click.option('--stranded_bigwigs', is_flag=True,
              help='Expect stranded, paired bigwigs with plus/forward '
              'strand first')
@click.option('--stranded_bed', is_flag=True, help='Expect stranded bed')
@click.option('--chrom_sizes', type=str, help='If specified, considers a '
              'chrom_sizes file when generating windows. If a window extends '
              'off a chromosome, the feature entry is skipped.')
def cli(bigwigs, feature_bed, output_matrix, anchor, bin_start, bin_number,
        bin_size, opposite_strand_fn, stranded_bigwigs, stranded_bed,
        chrom_sizes, temp_dir):
    """
    Generate matrices for stranded or unstranded bigWig matrices
    """
    BedMatrix(bigwigs, feature_bed, output_matrix, anchor, bin_start,
              bin_number, bin_size, opposite_strand_fn, stranded_bigwigs,
              stranded_bed, chrom_sizes, temp_dir)

if __name__ == '__main__':
    cli()
