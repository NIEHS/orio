import pymysql
import csv


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

    dbs = []
    for entry in cursor:
        dbs.append(entry[0])

    return dbs


def get_available_assemblies():
    assemblies = []

    dbs = get_databases()
    for db in dbs:
        cursor = get_UCSC_cursor()
        cursor.execute('SHOW TABLES IN {}'.format(db))
        tables = cursor.fetchall()
        for table in tables:
            if 'chromInfo' in table:
                assemblies.append(db)

    for assembly in ['hg19', 'mm9'][:-1]:
        assemblies.insert(0, assemblies.pop(assemblies.index(assembly)))

    return assemblies


def download_chromosome_sizes(path):
    dbs = get_databases()
    for db in dbs:
        cursor = get_UCSC_cursor()
        try:
            cursor.execute('SELECT chrom,size FROM {}.chromInfo'.format(db))
        except pymysql.InternalError as err:
            print("Error: {}".format(err))
        except pymysql.ProgrammingError as err:
            print("Error: {}".format(err))
        else:
            query_results = cursor.fetchall()

            fp = open(path + db + '.chromSizes.txt', 'w')
            chrom_sizes_file = csv.writer(fp, delimiter='\t')
            chrom_sizes_file.writerows(query_results)
            fp.close()
