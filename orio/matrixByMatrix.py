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


class MatrixByMatrix():

    def __init__(self, matrix_list, window_start,
                 bin_number, bin_size, sort_vector):

        self.matrix_list = matrix_list
        self.window_start = window_start
        self.bin_number = bin_number
        self.bin_size = bin_size
        self.sort_vector = sort_vector

        assert isinstance(window_start, int)
        assert isinstance(bin_number, int)
        assert isinstance(bin_size, int)

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
                            raise ValueError('headers not consistent across \
                                matrices')

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
                                raise ValueError('Row names do not match across \
                                    matrices')
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

        def getFCVectorData(vector_matrix, row_names, col_names, dsc_rep_data):
            fc_vectors = dict({'vectors': {}})

            if dsc_rep_data:
                fc_vectors['col_names'] = dsc_rep_data['col_names']

                for i, vector in enumerate(vector_matrix):
                    reorder_vector = []
                    for col_name in fc_vectors['col_names']:
                        reorder_vector.append(
                            vector_matrix[i][col_names.index(col_name)])
                    fc_vectors['vectors'][row_names[i]] = reorder_vector
            else:
                fc_vectors['col_names'] = col_names

                for i, vector in enumerate(vector_matrix):
                    fc_vectors['vectors'][row_names[i]] = vector

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

        def getFCCentroidData(kmeans_results, col_names, dsc_rep_data):
            fc_centroids = dict()

            for k in kmeans_results:
                fc_centroids[k] = dict()

                for cluster, centroid in \
                        enumerate(kmeans_results[k]['centroids']):

                    if dsc_rep_data is None:
                        fc_centroids[k][cluster+1] = centroid
                    else:
                        reorder_centroid = []
                        for col_name in dsc_rep_data['col_names']:
                            index = col_names.index(col_name)
                            reorder_centroid.append(centroid[index])
                        fc_centroids[k][cluster+1] = reorder_centroid

            return fc_centroids

        vector_matrix, row_names, col_names = \
            createVectorMatrix(self.matrix_files, self.matrix_names)
        kmeans_results = performKMeansClustering(vector_matrix)
        vector_matrix, kmeans_results = normalizeKMeans(
            vector_matrix, kmeans_results)
        self.fc_vectors = getFCVectorData(
            vector_matrix, row_names, col_names, self.dsc_rep_data)
        self.fc_clusters = getFCClusterData(
            kmeans_results, row_names)
        self.fc_centroids = getFCCentroidData(
            kmeans_results, col_names, self.dsc_rep_data)

    def minimizeMatrixOutputs(self):
        def reduceRowData(data):
            for i, row in enumerate(data['rows']):
                for j, value in enumerate(row['row_data']):
                    data['rows'][i]['row_data'][j] = \
                        '%.2f' % round(value, 2)

        def reduceCoordData(data):
            for i, coord in enumerate(data):
                for j, value in enumerate(coord):
                    data[i][j] = '%.2f' % round(value, 2)

        def reduceVectorData(data):
            for vector in data['vectors']:
                for i, value in enumerate(data['vectors'][vector]):
                    data['vectors'][vector][i] = \
                        '%.2f' % round(value, 2)

        def reduceCentroidData(data):
            for k in data:
                for cluster in data[k]:
                    for i, value in enumerate(data[k][cluster]):
                        data[k][cluster][i] = '%.2f' % round(value, 2)

        reduceRowData(self.dsc_full_data)
        reduceRowData(self.dsc_rep_data)

        reduceCoordData(self.dsc_dendrogram['icoord'])
        reduceCoordData(self.dsc_dendrogram['dcoord'])

        reduceVectorData(self.fc_vectors)
        reduceCentroidData(self.fc_centroids)

    def execute(self):
        self.readMatrixFiles()
        self.performDataSetClustering()
        self.performFeatureClustering()
        self.minimizeMatrixOutputs()

    def writeJson(self, fn):
        output_dict = {
            'dsc_full_data': self.dsc_full_data,
            'dsc_rep_data': self.dsc_rep_data,
            'dsc_dendrogram': self.dsc_dendrogram,
            'fc_vectors': self.fc_vectors,
            'fc_clusters': self.fc_clusters,
            'fc_centroids': self.fc_centroids,
        }
        with open(fn, 'w') as f:
            json.dump(output_dict, f, separators=(",", ": "))


@click.command()
@click.argument('matrix_list_fn', type=str)
@click.argument('window_start', type=int)
@click.argument('bin_number', type=int)
@click.argument('bin_size', type=int)
@click.argument('output_json', type=str)
@click.option('--sort_vector', nargs=1, type=str,
              help="Sort vector for correlative analysis")
def cli(matrix_list_fn, window_start, bin_number,
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
        matrix_list, window_start, bin_number, bin_size, sort_vector
    )
    mm.writeJson(output_json)


if __name__ == '__main__':
    cli()
