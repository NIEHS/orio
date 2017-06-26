import os

from .base import Validator
from .feature_list import FeatureListValidator


class SortVectorValidator(Validator):

    def __init__(self, feature_bed, sort_vector):

        super().__init__()

        assert os.path.exists(feature_bed)
        assert os.path.exists(sort_vector)

        self.feature_bed_fn = feature_bed
        self.sort_vector_fn = sort_vector

    def is_float(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    def check_if_parseable(self):
        try:
            with open(self.sort_vector_fn) as f:
                f.readlines()
        except UnicodeDecodeError:
            self.add_error('Invalid file format. Should be a tab-delimited '
                           'text format with two columns.')

    def check_columns(self):
        column_check = True
        float_check = True
        with open(self.sort_vector_fn) as f:
            for line in f:
                if len(line.split('\t')) != 2:
                    column_check = False
                if not self.is_float(line.strip().split('\t')[1]):
                    float_check = False

        if not column_check:
            self.add_error('Every row does not contain two columns')

        if not float_check:
            self.add_error('Column two contains a non-float value')

    def check_ids(self):
        with open(self.sort_vector_fn, 'r') as f:
            svs = f.read().splitlines()

        fbs = []
        with open(self.feature_bed_fn, 'r') as f:
            for line in f:
                if not FeatureListValidator.checkHeader(line):
                    fbs.append(line)

        if len(fbs) != len(svs):
            self.add_error('Feature list length ({}) does not equal sort vector length'.format(
                len(fbs), len(svs)))
            return

        for sv, fb in zip(svs, fbs):
            sv_id = sv.split('\t')[0]
            fl_id = fb.split('\t')[3]
            if sv_id != fl_id:
                self.add_error('Feature list ID ({}) does not match sort vector ID ({})'.format(
                    fl_id, sv_id))

    def validate(self):
        self.check_if_parseable()
        if self.validation_errors:
            return
        self.check_columns()
        self.check_ids()
