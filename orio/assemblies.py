import csv
import pymysql
import os
import re


def get_UCSC_cursor():
    cnx = pymysql.connect(
        user='genome',
        host='genome-mysql.cse.ucsc.edu',
    )
    cursor = cnx.cursor()
    return cursor


def get_databases():
    cursor = get_UCSC_cursor()
    cursor.execute('SHOW DATABASES')
    return [entry[0] for entry in cursor]


def get_assemblies():
    def bring_to_front(_list, value):
        _list.insert(0, _list.pop(_list.index(value)))

    cursor = get_UCSC_cursor()
    cursor.execute('SELECT table_schema from INFORMATION_SCHEMA.TABLES '
                   'where table_name="chromInfo";')
    assemblies = [d[0] for d in cursor.fetchall() if 'Patch' not in d[0]]

    for org in reversed(['hg', 'mm', 'dm', 'ce', 'sacCer']):
        org_list = [assembly for assembly in assemblies
                    if org == re.split('\d', assembly)[0]]
        for assembly in sorted(
            org_list, key=lambda x: int(re.split('\D+', x)[1])
        ):
            bring_to_front(assemblies, assembly)

    for assembly in reversed(['hg19', 'mm9']):
        bring_to_front(assemblies, assembly)

    return assemblies


def download_chromosome_sizes(path, overwrite=False):
    dbs = get_assemblies()
    for db in dbs:
        cursor = get_UCSC_cursor()
        try:
            cursor.execute('SELECT chrom, size FROM {}.chromInfo'.format(db))
        except (pymysql.InternalError, pymysql.ProgrammingError) as err:
            print("Error: {}".format(err))
        else:
            query_results = cursor.fetchall()
            fn = os.path.join(path, '{}.chromSizes.txt'.format(db))
            if overwrite or not os.path.isfile(fn):
                with open(fn, 'w') as f:
                    chrom_sizes_file = csv.writer(f, delimiter='\t')
                    chrom_sizes_file.writerows(query_results)
