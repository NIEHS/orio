#!/usr/bin/env python

import click
import numpy
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform, pdist
from scipy.cluster.vq import kmeans2, whiten
from math import floor
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

    def createSortOrders(self):
        self.sort_orders = {}

        for i, fn in enumerate(self.matrix_files):
            matrix_entries = dict()

            with open(fn) as f:
                next(f)
                for j, line in enumerate(f):
                    matrix_entries[j] = {
                        "feature_id": line.split()[0],
                        "row_sum": sum(map(float, line.strip().split()[1:]))
                    }

            order = sorted(
                matrix_entries,
                key=lambda x: (matrix_entries[x]["row_sum"]),
                reverse=True
            )

            self.sort_orders[self.matrix_ids[i]] = order

    def createVectorList(self):
        vector_list = []

        for fn in self.matrix_files:

            matrix = numpy.loadtxt(
                fn, skiprows=1, usecols=tuple(range(1, self.bin_number+1)),
                unpack=True
            )

            vector = numpy.sum(matrix, axis=0)
            vector_list.append(vector)

        return vector_list

    def createCorrelationMatrix(self):
        self.correlation_matrix = []

        if self.sort_vector is not None:
            sort_vector = self.readInSortVector(self.sort_vector)
            for matrix in self.matrix_files:
                self.correlation_matrix.append(
                    self.findVectorMatrixCorr(
                        sort_vector, matrix, self.bin_number
                    )
                )
            self.sort_vector = sort_vector
        else:
            vector_list = self.createVectorList()
            for i, vl1 in enumerate(vector_list):
                corrs = []
                for j, vl2 in enumerate(vector_list):
                    if i == j:
                        corr = 1.
                    else:
                        try:
                            corr = self.correlation_matrix[j][i]
                        except IndexError:
                            corr = stats.spearmanr(vl1, vl2)[0]
                    corrs.append(corr)
                self.correlation_matrix.append(corrs)

    def createDistanceMatrix(self):
        self.distance_matrix = []

        if self.sort_vector is not None:
            for cm1 in self.correlation_matrix:
                dists = []
                for cm2 in self.correlation_matrix:
                    dist = pdist(numpy.array([cm1, cm2]), "euclidean")[0]
                    dists.append(dist)
                self.distance_matrix.append(dists)
        else:
            for cm1 in self.correlation_matrix:
                dists = []
                for cm2 in cm1:
                    dists.append(1. - cm2)
                self.distance_matrix.append(dists)

    def findMedoids(self):
        self.cluster_medoids = []
        for cluster in self.cluster_members:
            if len(cluster) == 1:
                medoid = cluster[0]
            else:
                total_distances = []
                for i, clust1 in enumerate(cluster):
                    dist = 0.
                    index_1 = self.matrix_ids.index(clust1)
                    for j, clust2 in enumerate(cluster):
                        if i != j:
                            index_2 = self.matrix_ids.index(clust2)
                            dist += self.distance_matrix[index_1][index_2]
                    total_distances.append(dist)
                medoid = cluster[total_distances.index(min(total_distances))]
            self.cluster_medoids.append(medoid)

    def performClustering(self):
        self.createDistanceMatrix()

        distance_array = numpy.array(self.distance_matrix)
        lnk = linkage(squareform(distance_array), method="average")

        full_dg = dendrogram(lnk)
        truncated_dg = dendrogram(lnk, p=50, truncate_mode="lastp")

        self.dendrogram = truncated_dg
        self.cluster_members = []

        index = 0
        for entry in truncated_dg['ivl']:
            members = []
            if '(' in entry:
                count = int(entry.split("(")[1].split(")")[0])
                for i in range(count):
                    members.append(self.matrix_ids[int(full_dg['ivl'][index])])
                    index += 1
            else:
                members.append(self.matrix_ids[int(full_dg['ivl'][index])])
                index += 1
            self.cluster_members.append(members)
        self.findMedoids()

    def maxAbs(self, input_list):
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

    def findClusterCorrelationValues(self):

        self.cluster_correlation_values = []
        self.max_correlation_values = []
        self.med_correlation_values = []
        self.max_abs_correlation_values = []

        # define cluster_correlation_values
        for cm1 in self.cluster_members:
            ccvs = []
            if self.sort_vector is not None:
                for cm2 in cm1:
                    ccvs2 = []
                    index = self.matrix_ids.index(cm2)
                    for k, cm3 in enumerate(self.correlation_matrix[index]):
                        ccvs2.append(cm3)
                    ccvs.append(ccvs2)
            else:
                for cm2 in self.cluster_members:
                    ccvs2 = []
                    for cm3 in cm1:
                        for cm4 in cm2:
                            index_1 = self.matrix_ids.index(cm3)
                            index_2 = self.matrix_ids.index(cm4)
                            ccvs2.append(
                                self.correlation_matrix[index_1][index_2])
                    ccvs.append(ccvs2)
            self.cluster_correlation_values.append(ccvs)
        # get summary statistics
        if self.sort_vector is not None:
            for i, ccv1 in enumerate(self.cluster_correlation_values):
                sums = []
                for j, ccv2 in enumerate(ccv1):
                    sums.append(sum(ccv2))
                self.max_correlation_values.append(
                    ccv1[sums.index(max(sums))]
                )
                self.med_correlation_values.append(
                    ccv1[sums.index(sorted(sums)[int(floor((len(sums)-1)/2))])]
                )
                self.max_abs_correlation_values.append(
                    ccv1[sums.index(self.maxAbs(sums))]
                )
        else:
            for i, ccv1 in enumerate(self.cluster_correlation_values):

                max_cv = []
                med_cv = []
                max_abs_cv = []

                for j, ccv2 in enumerate(self.cluster_correlation_values[i]):
                    ccvs = self.cluster_correlation_values[i][j]
                    max_cv.append(max(ccvs))
                    med_cv.append(numpy.median(ccvs))
                    max_abs_cv.append(self.maxAbs(ccvs))

                self.max_correlation_values.append(max_cv)
                self.med_correlation_values.append(med_cv)
                self.max_abs_correlation_values.append(max_abs_cv)

    def getOutputDict(self):
        # make values more portable for output JSON
        for key in self.kmeans_results:
            for i in range(len(self.kmeans_results[key]['centroids'])):
                for j, val in enumerate(
                        self.kmeans_results[key]['centroids'][i]):
                    self.kmeans_results[key]['centroids'][i][j] = \
                        '%.2f' % round(val, 2)
        for entry in self.vector_matrix:
            for i, val in enumerate(entry):
                entry[i] = '%.2f' % round(val, 2)
        # Return an output dict of the analysis results
        if self.sort_vector is not None:
            self.sort_vector = self.sort_vector.tolist()
        return dict(
            bin_parameters={
                "window_start": self.window_start,
                "bin_number": self.bin_number,
                "bin_size": self.bin_size,
            },
            matrix_ids=self.matrix_ids,
            matrix_names=self.matrix_names,
            matrix_files=self.matrix_files,
            correlation_matrix=self.correlation_matrix,
            dendrogram=self.dendrogram,
            cluster_members=self.cluster_members,
            cluster_correlation_values=self.cluster_correlation_values,
            max_cluster_correlation_values=self.max_correlation_values,
            med_cluster_correlation_values=self.med_correlation_values,
            max_abs_correlation_values=self.max_abs_correlation_values,
            sort_orders=self.sort_orders,
            cluster_medoids=self.cluster_medoids,
            sort_vector=self.sort_vector,
            feature_clusters=self.kmeans_results,
            feature_vectors=self.vector_matrix,
            feature_columns=self.vector_columns,
            feature_names=self.feature_names,
            feature_cluster_members=self.feature_cluster_members,
        )

    def writeJson(self, fn):
        output_dict = self.getOutputDict()
        with open(fn, 'w') as f:
            json.dump(output_dict, f, separators=(",", ": "))

    def readInSortVector(self, input_file):
        """
        Loads sort vector using numpy
        """
        return numpy.loadtxt(input_file, usecols=[1])

    def findVectorMatrixCorr(self, vector, fn, bin_number):
        """
        Find correlation of each bin with vector in matrix
        """

        lst = []

        input_matrix = numpy.loadtxt(
            fn, skiprows=1,
            usecols=tuple(range(1, bin_number+1)),
            unpack=True)

        for i in range(len(input_matrix)):
            lst.append(stats.spearmanr(vector, input_matrix[i])[0])

        return lst

    def createFeatureMatrix(self):
        self.vector_matrix = None
        headers = None
        self.feature_names = []
        self.vector_columns = []
        for matrix in self.matrix_list:
            matrix_fn = matrix[2]
            self.vector_columns.append(matrix[1])
            with open(matrix_fn) as f:
                # DEAL WITH headersS
                # IF EMPTY, POPULATE headersS
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
                if not self.vector_matrix:
                    self.vector_matrix = []
                    for i, entry in enumerate(matrix_temp):
                        row_name = entry[0]
                        row_values = numpy.array(entry[1:]).astype(float)
                        self.feature_names.append(row_name)
                        self.vector_matrix.append([numpy.sum(row_values)])
                else:
                    for i, entry in enumerate(matrix_temp):
                        row_name = entry[0]
                        row_values = numpy.array(entry[1:]).astype(float)
                        if row_name != self.feature_names[i]:
                            raise ValueError('Row names do not match across \
                                matrices')
                        self.vector_matrix[i].append(numpy.sum(row_values))

    def performFeatureClustering(self):
        self.kmeans_results = dict()
        whitened = whiten(self.vector_matrix)
        std_devs = numpy.std(self.vector_matrix, axis=0)
        for k in range(2, 11):
            centroids, labels = kmeans2(whitened, k, minit='points')
            self.kmeans_results[k] = {
                'centroids': centroids.tolist(),
                'labels': labels.tolist()
                }
            for i, centroid in enumerate(self.kmeans_results[k]['centroids']):
                for j, val in enumerate(centroid):
                    self.kmeans_results[k]['centroids'][i][j] = \
                        val * std_devs[j]

    def normalizeFeatureMatrix(self):
        # Normalize feature vectors
        vnorm = numpy.percentile(self.vector_matrix, 75, axis=0)

        for i, vector in enumerate(self.vector_matrix):
            for j, val in enumerate(vector):
                self.vector_matrix[i][j] = \
                    val / vnorm[j]

        for k in self.kmeans_results:
            for i, centroid in enumerate(self.kmeans_results[k]['centroids']):
                for j, val in enumerate(centroid):
                    self.kmeans_results[k]['centroids'][i][j] = \
                        val / vnorm[j]

    def getClusterMembers(self):
        self.feature_cluster_members = dict()
        for k in self.kmeans_results:
            self.feature_cluster_members[k] = dict()
            for i in range(1, k+1):
                self.feature_cluster_members[k][i] = []
            for i, cluster in enumerate(self.kmeans_results[k]['labels']):
                self.feature_cluster_members[k][int(cluster)+1].append(
                    self.feature_names[i])

    def reorderFeatureMatrixByDendrogram(self):
        order = []
        for entry in self.cluster_medoids:
            order.append(self.matrix_ids.index(entry))
        for i, row in enumerate(self.vector_matrix):
            temp = []
            for index in order:
                temp.append(row[index])
            self.vector_matrix[i] = temp

    def execute(self):
        self.readMatrixFiles()
        self.createSortOrders()
        self.createCorrelationMatrix()
        self.performClustering()
        self.findClusterCorrelationValues()

        self.createFeatureMatrix()
        self.reorderFeatureMatrixByDendrogram()
        self.performFeatureClustering()
        self.normalizeFeatureMatrix()
        self.getClusterMembers()


@click.command()
@click.argument('matrix_list_fn', type=str)
@click.argument('window_start', type=int)
@click.argument('bin_number', type=int)
@click.argument('bin_size', type=int)
@click.argument('output_json', type=str)
@click.option('--sort_vector', nargs=1, type=str, help="Sort vector for correlative analysis")
def cli(matrix_list_fn, window_start, bin_number,
        bin_size, output_json, sort_vector):
    """
    Considering matrix files specified by a list, performs cross-matrix
    correlative analysis.

    \b
    Arguments:
    - matrix_list_fn:   List of matrix files to be considered in analysis. Each
                        row in the list corresponds to a matrix to be considered
                        in the analysis. The list contains three columns:
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

    mm = MatrixByMatrix(matrix_list, window_start, bin_number, bin_size, sort_vector)
    mm.writeJson(output_json)


if __name__ == '__main__':
    cli()
