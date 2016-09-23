import csv
import pymysql
import os

from . import utils


def get_ucsc_cursor(cursorclass=pymysql.cursors.Cursor):
    cnx = pymysql.connect(
        user='genome',
        host='genome-mysql.cse.ucsc.edu',
        cursorclass=cursorclass,
    )
    cursor = cnx.cursor()
    return cursor


def get_databases():
    cursor = get_ucsc_cursor()
    cursor.execute('SHOW DATABASES')
    return [entry[0] for entry in cursor]


def get_assemblies():
    cursor = get_ucsc_cursor()
    cursor.execute('SELECT table_schema from INFORMATION_SCHEMA.TABLES where table_name="chromInfo";')  # noqa
    return [d[0] for d in cursor.fetchall()]


def download_annotations():
    assemblies = [
        'hg19',
        'mm9',
        'hg38',
        'mm10',
        'dm6',
    ]
    cursor = get_ucsc_cursor(cursorclass=pymysql.cursors.DictCursor)
    for assembly in assemblies:
        query = """
        SELECT {0}.refGene.name, {0}.refGene.chrom, {0}.refGene.exonStarts,
               {0}.refGene.exonEnds, {0}.refGene.strand, {0}.refGene.name2,
               hgFixed.refLink.locusLinkId
        FROM {0}.refGene LEFT OUTER JOIN
            hgFixed.refLink ON {0}.refGene.name=hgFixed.refLink.mrnaAcc
        """.format(assembly)
        cursor.execute(query)
        transcripts = [entry for entry in cursor]

        sorted_transcripts = sorted(transcripts, key=lambda x: (
            x['chrom'],
            int(x['exonStarts'].decode('UTF-8').split(',')[0]),
            x['name'].split('_')[0],
            x['name'].split('_')[1]
        ))

        path = utils.get_data_path()
        fn = os.path.join(path, '{}.gtf'.format(assembly))
        with open(fn, 'w') as f:
            for entry in sorted_transcripts:
                starts = entry['exonStarts'].decode('UTF-8').split(',')[:-1]
                ends = entry['exonEnds'].decode('UTF-8').split(',')[:-1]
                for start, end in zip(starts, ends):
                    f.write('\t'.join([
                        entry['chrom'],
                        assembly + '.refGene',
                        'exon',
                        start,
                        end,
                        '.',
                        entry['strand'],
                        '.',
                        'gene_id "{}"; transcript_id "{}";'.format(entry['name2'], entry['name'],),
                    ]) + '\n')

    cursor.close()


def download_chromosome_sizes(path):
    dbs = get_assemblies()
    for db in dbs:
        cursor = get_ucsc_cursor()
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
