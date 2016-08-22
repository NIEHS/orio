#!/usr/bin/env python

import click
import numpy
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform, pdist
from scipy.cluster.vq import kmeans2, whiten
import copy
import os
import json
import bisect
from collections import defaultdict


def readGTF(annotation_file):
    transcripts = defaultdict(dict)

    with open(annotation_file) as f:
        for line in f:
            chromosome, source, feature, start, end, score, strand, frame,\
                attributes = line.strip().split('\t')

            keys = []
            values = []
            gtf_fields = dict()
            for entry in attributes.split(';')[:-1]:
                keys.append(entry.split('\"')[0].strip())
                values.append(entry.split('\"')[1].strip())
            for key, value in zip(keys, values):
                gtf_fields[key] = value

            tr_id = gtf_fields.pop('transcript_id')
            gene_id = gtf_fields.pop('gene_id')

            if feature == 'exon':
                if tr_id not in transcripts[chromosome]:
                    transcripts[chromosome][tr_id] = {
                        'strand': strand,
                        'exons': [],
                        'gene_id': gene_id,
                        }
                transcripts[chromosome][tr_id]['exons'].append(
                    [int(start), int(end)]
                    )

    for chromosome in transcripts:
        for transcript in transcripts[chromosome].values():
            transcript['exons'].sort(key=lambda x: x[0])
            transcript['start'] = int(transcript['exons'][0][0])
            transcript['end'] = int(transcript['exons'][-1][1])

    return transcripts


def checkHeader(line):
    # Check to see if line is header
    if line == "\n":
        return True
    elif line[0] == "#":
        return True
    elif line.split()[0].lower() in ("track", "browser"):
        return True
    else:
        return False


def readBED(bed_fn):
    bed_list = []
    with open(bed_fn) as bed:
        for line in bed:
            if not checkHeader(line):
                chromosome, start, end = line.strip().split()[0:3]
                bed_list.append({
                    'chromosome': chromosome,
                    'start': int(start) + 1,
                    'end': int(end),
                })
    return bed_list


def associateGenes(features, transcripts):
    gene_list = []

    start_values = dict()
    start_keys = dict()

    end_values = dict()
    end_keys = dict()

    for chromosome in transcripts:
        sort = \
            sorted(transcripts[chromosome].items(),
                   key=lambda x: x[1]['start'])
        start_values[chromosome] = \
            [r[1]['start'] for r in sort]
        start_keys[chromosome] = \
            [r[0] for r in sort]

        sort = \
            sorted(transcripts[chromosome].items(),
                   key=lambda x: x[1]['end'])
        end_values[chromosome] = \
            [r[1]['end'] for r in sort]
        end_keys[chromosome] = \
            [r[0] for r in sort]

    for feature in features:
        chromosome = feature['chromosome']

        left_bisect = bisect.bisect_left(
            end_values[chromosome],
            feature['start']
        )
        right_bisect = bisect.bisect_right(
            start_values[chromosome],
            feature['end']
        )
        left_index = max(left_bisect-1, 0)
        right_index = min(right_bisect+1, len(start_values[chromosome]))

        shortest_dist = float('Inf')
        closest_gene = set()

        for transcript in \
                set(start_keys[chromosome][left_index:]) & \
                set(end_keys[chromosome][0:right_index]):

            dist = numpy.amax([
                0,
                feature['start'] - transcripts[chromosome][transcript]['end'],
                transcripts[chromosome][transcript]['start'] - feature['end'],
            ])
            if dist < shortest_dist:
                shortest_dist = dist
                closest_gene = \
                    {transcripts[chromosome][transcript]['gene_id']}
            elif dist == shortest_dist:
                closest_gene.add(
                    transcripts[chromosome][transcript]['gene_id']
                )

        gene_list.append(','.join(list(closest_gene)))

    return gene_list


class MatrixByMatrix():

    def __init__(self, feature_bed, matrix_list, annotation, window_start,
                 bin_number, bin_size, sort_vector):

        self.feature_bed = feature_bed
        self.annotation = annotation
        self.matrix_list = matrix_list
        self.window_start = window_start
        self.bin_number = bin_number
        self.bin_size = bin_size
        self.sort_vector = sort_vector

        assert isinstance(window_start, int)
        assert isinstance(bin_number, int)
        assert isinstance(bin_size, int)

        assert os.path.exists(self.annotation)
        assert os.path.exists(self.feature_bed)
        if sort_vector:
            assert os.path.exists(self.sort_vector)

        self.execute()

    def readMatrixFiles(self):
        self.matrix_ids, \
            self.matrix_names, \
            self.matrix_files,  = zip(*self.matrix_list)

    def performDataSetClustering(self):

        def readInSortVector(input_fn):
            return numpy.loadtxt(input_fn, usecols=[1])

        def createVectorList(file_list, bin_number):
            vector_list = []

            for fn in file_list:
                matrix = numpy.loadtxt(
                    fn, skiprows=1,
                    usecols=tuple(range(1, bin_number+1)),
                    unpack=True
                )
                vector = numpy.sum(matrix, axis=0)
                vector_list.append(vector)

            return vector_list

        def findVectorMatrixCorr(vector, fn, bin_number):
            lst = []
            input_matrix = numpy.loadtxt(
                fn, skiprows=1,
                usecols=tuple(range(1, bin_number+1)),
                unpack=True)
            for i in range(len(input_matrix)):
                lst.append(stats.spearmanr(vector, input_matrix[i])[0])
            return lst

        def createCorrelationMatrix(matrix_files, sort_vector, bin_number):

            corr_matrix = []

            if sort_vector is not None:
                sort_vector = readInSortVector(sort_vector)
                for matrix in matrix_files:
                    corr_matrix.append(
                        findVectorMatrixCorr(
                            sort_vector, matrix, bin_number
                        )
                    )
                sort_vector = sort_vector
            else:
                vector_list = createVectorList(matrix_files, bin_number)
                for i, vl1 in enumerate(vector_list):
                    corrs = []
                    for j, vl2 in enumerate(vector_list):
                        if i == j:
                            corr = 1.
                        else:
                            try:
                                corr = corr_matrix[j][i]
                            except IndexError:
                                corr = stats.spearmanr(vl1, vl2)[0]
                        corrs.append(corr)
                    corr_matrix.append(corrs)
            return corr_matrix

        def createDistanceMatrix(corr_matrix, sort_vector):
            dist_matrix = []

            if sort_vector is not None:
                for cm1 in corr_matrix:
                    dists = []
                    for cm2 in corr_matrix:
                        dist = \
                            pdist(numpy.array([cm1, cm2]), "euclidean")[0]
                        dists.append(dist)
                    dist_matrix.append(dists)
            else:
                for cm1 in corr_matrix:
                    dists = []
                    for cm2 in cm1:
                        dists.append(1. - cm2)
                    dist_matrix.append(dists)

            return dist_matrix

        def getClusterMembers(truncated_dg, full_dg, matrix_ids):
            cluster_members = []

            index = 0
            for entry in truncated_dg['ivl']:
                members = []
                if '(' in entry:
                    count = int(entry.split("(")[1].split(")")[0])
                    for i in range(count):
                        members.append(int(full_dg['ivl'][index]))
                        index += 1
                else:
                    members.append(int(full_dg['ivl'][index]))
                    index += 1
                cluster_members.append(members)

            return cluster_members

        def getClusterMedoids(cluster_members, matrix_ids, dist_matrix):
            cluster_medoids = []

            for cluster in cluster_members:
                if len(cluster) == 1:
                    medoid = cluster[0]
                else:
                    total_distances = []
                    for i in cluster:
                        dist = 0.
                        for j in cluster:
                            if i != j:
                                dist += dist_matrix[i][j]
                        total_distances.append(dist)
                    medoid = \
                        cluster[total_distances.index(min(total_distances))]
                cluster_medoids.append(medoid)

            return cluster_medoids

        def getBinNames(window_start, bin_number, bin_size):
            bin_names = []

            for i in range(bin_number):
                bin_names.append(str(window_start + i * bin_size) +
                                 ':' +
                                 str((window_start + (i+1) * bin_size)-1))

            return bin_names

        def getFullCorrData(cluster_members, matrix_ids, matrix_names,
                            corr_matrix, sort_vector, bin_names):
            full_corr = dict()
            matrix_order = []

            for cluster in cluster_members:
                for id_num in cluster:
                    matrix_order.append(id_num)

            full_corr.update({'col_names': [], 'rows': []})
            if sort_vector is None:
                full_corr['col_ids'] = []

            for i in matrix_order:
                if sort_vector is None:
                    full_corr['col_ids'].append(matrix_ids[i])
                    full_corr['col_names'].append(matrix_names[i])
                else:
                    full_corr['col_names'] = bin_names

                if sort_vector is None:
                    row_data = []
                    for j in matrix_order:
                        row_data.append(corr_matrix[i][j])
                else:
                    row_data = corr_matrix[i]

                full_corr['rows'].append({
                    'row_id': matrix_ids[i],
                    'row_name': matrix_names[i],
                    'row_data': row_data,
                })

            return full_corr

        def maxAbs(input_list):
            abs_list = []
            for entry in input_list:
                abs_list.append(abs(entry))
            max_index = None
            max_value = None
            for i in range(len(abs_list)):
                if max_index is None:
                    max_index = i
                    max_value = abs_list[i]
                elif abs_list[i] > max_value:
                    max_index = i
                    max_value = abs_list[i]
                elif abs_list[i] == max_value:
                    if input_list[i] > input_list[max_index]:
                        max_index = i
            return input_list[max_index]

        def getClusterRepRow(index, cluster_members, corr_matrix, sort_vector):
            if sort_vector:
                rows = []
                row_sums = []
                for i in cluster_members[index]:
                    rows.append(corr_matrix[i])
                    row_sums.append(sum(rows[-1]))
                return rows[row_sums.index(max(row_sums))]
            else:
                row_data = []
                for cluster in cluster_members:
                    row_values = []
                    for i in cluster:
                        for j in cluster_members[index]:
                            row_values.append(corr_matrix[i][j])
                    row_data.append(maxAbs(row_values))
                return row_data

        def getRepCorrData(cluster_medoids, cluster_members, matrix_ids,
                           matrix_names, corr_matrix, sort_vector, bin_names):
            rep_corr = dict()

            rep_corr.update({'col_names': [], 'rows': []})
            if sort_vector is None:
                rep_corr['col_ids'] = []

            for clust_num, index in enumerate(cluster_medoids):
                if sort_vector is None:
                    rep_corr['col_ids'].append(matrix_ids[index])
                    rep_corr['col_names'].append(matrix_names[index])
                else:
                    rep_corr['col_names'] = bin_names

                row_data = getClusterRepRow(
                    clust_num, cluster_members, corr_matrix, sort_vector)

                rep_corr['rows'].append({
                    'row_id': matrix_ids[index],
                    'row_name': matrix_names[index],
                    'row_data': row_data,
                })

            return rep_corr

        corr_matrix = createCorrelationMatrix(
            self.matrix_files, self.sort_vector, self.bin_number)
        dist_matrix = createDistanceMatrix(corr_matrix, self.sort_vector)

        dist_array = numpy.array(dist_matrix)
        lnk = linkage(squareform(dist_array), method="average")

        full_dg = dendrogram(lnk)
        truncated_dg = dendrogram(lnk, p=50, truncate_mode="lastp")

        self.dsc_dendrogram = truncated_dg

        cluster_members = getClusterMembers(
            truncated_dg, full_dg, self.matrix_ids)
        cluster_medoids = getClusterMedoids(
            cluster_members, self.matrix_ids, dist_matrix)

        bin_names = getBinNames(
            self.window_start, self.bin_number, self.bin_size)

        self.dsc_full_data = getFullCorrData(
            cluster_members, self.matrix_ids, self.matrix_names, corr_matrix,
            self.sort_vector, bin_names)
        self.dsc_rep_data = getRepCorrData(
            cluster_medoids, cluster_members, self.matrix_ids,
            self.matrix_names, corr_matrix, self.sort_vector, bin_names)

    def performFeatureClustering(self):
        def readMatrixOrder(dsc_rep_data):
            matrix_order = []

            for row in dsc_rep_data['rows']:
                matrix_order.append(row['row_name'])

            return matrix_order

        def createVectorMatrix(matrix_files, matrix_names):
            vector_matrix = None
            headers = None
            row_names = []
            col_names = []

            for matrix_fn, matrix_name in zip(matrix_files, matrix_names):
                col_names.append(matrix_name)
                with open(matrix_fn) as f:
                    # DEAL WITH HEADERS
                    # IF EMPTY, POPULATE HEADERS
                    if not headers:
                        headers = next(f).strip().split()
                    # ELSE, CHECK IF CONSISTENT
                    else:
                        if headers != next(f).strip().split():
                            raise ValueError('headers not consistent across matrices')  # noqa

                    # POPULATE TEMPORARY MATRIX
                    matrix_temp = []
                    for line in f:
                        matrix_temp.append(line.strip().split())

                    # ADD SUM TO VECTOR MATRIX
                    if not vector_matrix:
                        vector_matrix = []
                        for i, entry in enumerate(matrix_temp):
                            row_name = entry[0]
                            row_values = numpy.array(entry[1:]).astype(float)
                            row_names.append(row_name)
                            vector_matrix.append([numpy.sum(row_values)])
                    else:
                        for i, entry in enumerate(matrix_temp):
                            row_name = entry[0]
                            row_values = numpy.array(entry[1:]).astype(float)
                            if row_name != row_names[i]:
                                raise ValueError('Row names do not match across matrices')  # noqa
                            vector_matrix[i].append(numpy.sum(row_values))

            return vector_matrix, row_names, col_names

        def performKMeansClustering(vector_matrix):
            kmeans_results = dict()
            whitened = whiten(vector_matrix)
            std_devs = numpy.std(vector_matrix, axis=0)
            for k in range(2, 11):
                centroids, labels = kmeans2(whitened, k, minit='points')
                kmeans_results[k] = {
                    'centroids': centroids.tolist(),
                    'labels': labels.tolist()
                    }
                for i, centroid in enumerate(kmeans_results[k]['centroids']):
                    for j, val in enumerate(centroid):
                        kmeans_results[k]['centroids'][i][j] = \
                            val * std_devs[j]
            return kmeans_results

        def normalizeKMeans(vector_matrix, kmeans_results):
            vnorm = numpy.percentile(vector_matrix, 75, axis=0)
            norm_matrix = copy.copy(vector_matrix)
            norm_kmeans = copy.copy(kmeans_results)

            for i, vector in enumerate(vector_matrix):
                for j, val in enumerate(vector):
                    norm_matrix[i][j] = \
                        val / vnorm[j]

            for k in kmeans_results:
                for i, centroid in enumerate(kmeans_results[k]['centroids']):
                    for j, val in enumerate(centroid):
                        norm_kmeans[k]['centroids'][i][j] = \
                            val / vnorm[j]

            return norm_matrix, norm_kmeans

        def getFCVectorData(vector_matrix, row_names, col_names, matrix_order):
            fc_vectors = dict({'vectors': {}})

            # if dsc_rep_data and not sort_vector:
            fc_vectors['col_names'] = matrix_order

            for i, vector in enumerate(vector_matrix):
                reorder_vector = []
                for matrix in matrix_order:
                    reorder_vector.append(
                        vector_matrix[i][col_names.index(matrix)])
                fc_vectors['vectors'][row_names[i]] = reorder_vector

            return fc_vectors

        def getFCClusterData(kmeans_results, row_names):
            fc_clusters = dict()

            for k in kmeans_results:
                fc_clusters[k] = dict()
                for i in range(1, k+1):
                    fc_clusters[k].update({i: []})
                for i, cluster in enumerate(kmeans_results[k]['labels']):
                    fc_clusters[k][cluster+1].append(row_names[i])

            return fc_clusters

        def getFCCentroidData(kmeans_results, col_names, matrix_order):
            fc_centroids = dict()

            for k in kmeans_results:
                fc_centroids[k] = dict()

                for cluster, centroid in enumerate(kmeans_results[k]['centroids']):  # noqa

                    reorder_centroid = []
                    for matrix in matrix_order:
                        index = col_names.index(matrix)
                        reorder_centroid.append(centroid[index])
                    fc_centroids[k][cluster+1] = reorder_centroid

            return fc_centroids

        matrix_order = readMatrixOrder(self.dsc_rep_data)
        vector_matrix, row_names, col_names = \
            createVectorMatrix(self.matrix_files, self.matrix_names)
        kmeans_results = performKMeansClustering(vector_matrix)
        vector_matrix, kmeans_results = normalizeKMeans(
            vector_matrix, kmeans_results)

        self.fc_vectors = getFCVectorData(
            vector_matrix, row_names, col_names, matrix_order)
        self.fc_clusters = getFCClusterData(
            kmeans_results, row_names)
        self.fc_centroids = getFCCentroidData(
            kmeans_results, col_names, matrix_order)

    def readRowNames(self, matrix_fn):
        row_names = []
        with open(matrix_fn) as f:
            next(f)
            for line in f:
                name = line.strip().split()[0]
                row_names.append(name)
        return row_names

    def findClosestGene(self):
        features = readBED(self.feature_bed)
        genes = readGTF(self.annotation)

        feature_list = self.readRowNames(self.matrix_list[0][2])
        gene_list = associateGenes(features, genes)

        self.feature_to_gene = {k: v for k, v in zip(feature_list, gene_list)}

    def execute(self):
        self.readMatrixFiles()
        self.performDataSetClustering()
        self.performFeatureClustering()
        self.findClosestGene()

    def writeJson(self, fn):
        output_dict = {
            'dsc_full_data': self.dsc_full_data,
            'dsc_rep_data': self.dsc_rep_data,
            'dsc_dendrogram': self.dsc_dendrogram,
            'fc_vectors': self.fc_vectors,
            'fc_clusters': self.fc_clusters,
            'fc_centroids': self.fc_centroids,
            'feature_to_gene': self.feature_to_gene,
        }
        with open(fn, 'w') as f:
            json.dump(output_dict, f, separators=(",", ": "))


@click.command()
@click.argument('feature_bed', type=str)
@click.argument('matrix_list_fn', type=str)
@click.argument('annotation', type=str)
@click.argument('window_start', type=int)
@click.argument('bin_number', type=int)
@click.argument('bin_size', type=int)
@click.argument('output_json', type=str)
@click.option('--sort_vector', nargs=1, type=str,
              help="Sort vector for correlative analysis")
def cli(feature_bed, matrix_list_fn, annotation, window_start, bin_number,
        bin_size, output_json, sort_vector):
    """
    Considering matrix files specified by a list, performs cross-matrix
    correlative analysis.

    \b
    Arguments:
    - matrix_list_fn:   List of matrix files to be considered in analysis. Each
                        row in the list corresponds to a matrix to be
                        considered in the analysis. The list contains three
                        columns:
                            1) unique integer ID for matrix
                            2) unique name for each matrix
                            3) absolute path to matrix file
    - window_start:     The first position of the analysis window relative to
                        the features in the associated feature list.
    - bin_number:       Number of bins used
    - bin_size:         Size of each bin
    - output_json:      Filename of output JSON
    """

    assert os.path.exists(matrix_list_fn)
    with open(matrix_list_fn) as f:
        matrix_list = [
            line.strip().split()
            for line in f.readlines()
        ]

    mm = MatrixByMatrix(
        feature_bed, matrix_list, annotation, window_start, bin_number,
        bin_size, sort_vector
    )
    mm.writeJson(output_json)


if __name__ == '__main__':
    cli()
