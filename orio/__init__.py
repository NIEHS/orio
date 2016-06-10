from .downloadBinaries import binaries_exist, download_ucsc_tools

__version__ = '0.0.1'


if not binaries_exist():
    download_ucsc_tools()
