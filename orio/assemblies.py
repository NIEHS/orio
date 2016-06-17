import csv
import pymysql
import os


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
    cursor = get_UCSC_cursor()
    cursor.execute('SELECT table_schema from INFORMATION_SCHEMA.TABLES where table_name="chromInfo";')
    return [d[0] for d in cursor.fetchall()]


def download_chromosome_sizes(path):
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
            with open(fn, 'w') as f:
                chrom_sizes_file = csv.writer(f, delimiter='\t')
                chrom_sizes_file.writerows(query_results)
