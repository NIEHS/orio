import os
from abc import abstractmethod

from .. import utils


class Validator(object):

    def __init__(self):
        self.validation_errors = []
        self.validation_warnings = []

    @property
    def is_valid(self):
        return len(self.validation_errors) == 0

    @abstractmethod
    def validate(self):
        pass

    def add_error(self, txt):
        self.validation_errors.append(txt)

    def add_errors(self, lst):
        self.validation_errors.extend(lst)

    def add_warning(self, txt):
        self.validation_warnings.append(txt)

    def add_warnings(self, lst):
        self.validation_warnings.extend(lst)

    def display_errors(self):
        return '\n'.join(set(self.validation_errors))

    def display_warnings(self):
        return '\n'.join(set(self.validation_warnings))


def get_validate_files_path():
    path = os.path.join(utils.get_bin_path(), 'validateFiles')
    if not os.path.exists(path):
        raise IOError('validateFiles not found, expected {}'.format(path))
    return path


def get_chromosome_size_path(genome):
    path = os.path.join(utils.get_data_path(), genome + '.chromSizes.txt')
    if not os.path.exists(path):
        raise IOError('File not found: {0}'.format(path))
    return path


def get_annotation_file_path(genome):
    path = os.path.join(utils.get_data_path(), genome + '.gtf')
    if not os.path.exists(path):
        raise IOError('File not found: {0}'.format(path))
    return path
