from .downloadBinaries import binaries_exist, download_ucsc_tools  # noqa
import os

__version__ = '0.0.1'


def get_data_path():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))


def get_bin_path():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), 'bin'))
