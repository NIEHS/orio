import os
from abc import abstractmethod


class Validator(object):

    def __init__(self):
        self.validation_errors = []

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

    def display_errors(self):
        return '\n'.join(set(self.validation_errors))


def get_validateFiles_path():
    path = os.path.abspath(
        os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            os.path.pardir,
            'validateFiles'
        )
    )
    if not os.path.exists(path):
        raise IOError('validateFiles not found, expected {}'.format(path))
    return path


def get_chromosome_size_path(genome):
    root = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(root, 'data', genome + '.chromSizes.txt')
    if not os.path.exists(path):
        raise IOError('File not found: {0}'.format(path))
    return path
