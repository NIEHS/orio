import os

from .base import Validator


class FeatureListValidator(Validator):

    # Possible 'empty' names used in a BED file
    DUMMY_NAMES = ['.']

    def __init__(self, feature_list, chrom_sizes_file, stranded):

        super().__init__()

        assert os.path.exists(feature_list)
        assert os.path.exists(chrom_sizes_file)

        self.feature_list = feature_list
        self.chrom_sizes_file = chrom_sizes_file
        self.stranded = stranded

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
    def readChromSizes(sizes_file):
        sizes = dict()
        with open(sizes_file) as f:
            for line in f:
                chrom, size = line.strip().split()
                sizes[chrom] = int(size)
        return sizes

    def check_if_parseable(self):
        try:
            with open(self.feature_list) as f:
                f.readlines()
        except UnicodeDecodeError:
            self.add_error('Invalid file format. Should be in BED text-file format.')

    def set_number_columns(self):
        # Find number of columns in bed
        self.number_columns = None
        with open(self.feature_list) as f:
            for line in f:
                if not (self.checkHeader(line)):
                    if self.number_columns:
                        if self.number_columns != len(line.split()):
                            self.add_error('Inconsistent number of columns')
                    else:
                        self.number_columns = len(line.split())

    def check_unique_feature_names(self):
        # If BED file contains names (cols >= 4), make sure they are unique
        if self.number_columns < 4:
            return

        feature_names = set()
        with open(self.feature_list) as f:
            contains_names = False
            contains_dummy = False

            for line in f:
                if not (self.checkHeader(line)):
                    feature_name = line.strip().split()[3]
                    if feature_name in self.DUMMY_NAMES:
                        contains_dummy = True
                    else:
                        contains_names = True
                        if feature_name in feature_names:
                            self.add_error('Duplicate feature name: {}'.format(feature_name))  # noqa
                        else:
                            feature_names.add(feature_name)

            if contains_names and contains_dummy:
                self.add_error('Features inconsistently have names')

    def check_bed_entries(self):
        sizes = self.readChromSizes(self.chrom_sizes_file)
        with open(self.feature_list) as f:
            for line in f:
                if not self.checkHeader(line):
                    # Check for three columns
                    try:
                        chrom, start, end = line.strip().split()[0:3]
                    except(ValueError):
                        self.add_error('Too few columns: {}'.format(line.strip()))  # noqa
                    else:
                        # Check if chrom in sizes file
                        if chrom not in sizes:
                            self.add_error('Chromosome "{}" not in chromosome list'.format(chrom))  # noqa
                        # Check if entry is contained within chromosome
                        elif int(start) > sizes[chrom] or \
                                int(end) > sizes[chrom] or \
                                int(start) < 0 or \
                                int(end) < 1:
                            self.add_error('Entry not within chromosome: {}'.format(line.strip()))  # noqa

    def check_stranded(self):
        if self.stranded:
            with open(self.feature_list) as f:
                for line in f:
                    if not self.checkHeader(line):
                        try:
                            strand = line.strip().split()[5]
                        except(IndexError):
                            self.add_error('Missing strand column in feature list: {}'.format(line.strip()))  # noqa
                        else:
                            if strand != '+' and strand != '-':
                                self.add_error('Strand not \'+\' or \'-\': {}'.format(line.strip()))  # noqa

    def validate(self):
        self.check_if_parseable()
        if self.validation_errors:
            return
        self.set_number_columns()
        self.check_bed_entries()
        self.check_unique_feature_names()
        self.check_stranded()
