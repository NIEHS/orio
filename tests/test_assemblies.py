import pytest

from orio import utils
from orio import assemblies


def test_cursor():
    cursor = assemblies.get_UCSC_cursor()
    cursor.close()


def test_get_databases():
    dbs = assemblies.get_databases()
    assert 'hg19' in dbs
    assert 'mm9' in dbs
    assert len(dbs) > 200


def test_get_assemblies():
    dbs = assemblies.get_assemblies()
    assert 'hg19' in dbs
    assert 'mm9' in dbs
    assert len(dbs) >= 190


@pytest.mark.skip(reason="slow")
def test_download_chromosome_sizes():
    path = utils.get_data_path()
    assemblies.download_chromosome_sizes(path)


@pytest.mark.skip(reason="slow")
def test_download_annotations():
    assemblies.download_annotations()
