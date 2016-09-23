#!/usr/bin/env python

import click
import numpy
from scipy.cluster.vq import kmeans2, whiten
import os
import json


class ClusterFeatures():

    def __init__(self, matrix_list):
        self.matrix_list = matrix_list
        self.execute()

    def create_vector_matrix(self):
        headers = None
        vector_matrix = None
        row_names = []

        for entry in self.matrix_list:
            matrix_fn = entry[2]
            with open(matrix_fn) as f:

                # Populate headers if empty
                if headers is None:
                    headers = next(f).strip().split()
                # ...otherwise ensure headers are consistent
                else:
                    if headers != next(f).strip().split():
                        raise ValueError('Headers not consistent across matrices')

                # Create a temporary matrix
                matrix_temp = []
                for line in f:
                    matrix_temp.append(line.strip().split())

                # Add sum to vector matrix
                if vector_matrix is None:
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
                            raise ValueError('Row names do not match across matrices')
                        vector_matrix[i].append(numpy.sum(row_values))

        self.headers = headers
        self.row_names = row_names
        self.vector_matrix = vector_matrix

    def perform_clustering(self):
        whitened = whiten(self.vector_matrix)
        for i in range(2, 11):
            centroids, labels = kmeans2(whitened, i)
            self.kmeans_results[i] = dict(
                centroids=centroids.tolist(),
                labels=labels.tolist()
            )

    def write_json(self, fn):
        d = dict(
            kmeans_results=self.kmeans_results,
            vector_matrix=self.vector_matrix,
            bins=self.headers,
            row_names=self.row_names
        )
        with open(fn, 'w') as f:
            json.dump(d, f, separators=(",", ": "))

    def execute(self):
        self.kmeans_results = {}
        self.create_vector_matrix()
        self.perform_clustering()


@click.command()
@click.argument('matrix_list_fn', type=str)
@click.argument('output_json', type=str)
def cli(matrix_list_fn, output_json):
    r"""
    Cluster features (rows) by vectors derived from matrix list.

    \b
    Arguments:
    - matrix_list_fn:   List of matrix files to be considered in analysis. Each
                        row in the list corresponds to a matrix to be considered
                        in the analysis. The list contains three columns:
                            1) unique integer ID for matrix
                            2) unique name for each matrix
                            3) absolute path to matrix file
    - output_json:      Filename of output JSON
    """
    assert os.path.exists(matrix_list_fn)
    with open(matrix_list_fn) as f:
        matrix_list = [
            line.strip().split()
            for line in f.readlines()
        ]

    cf = ClusterFeatures(matrix_list)
    cf.write_json(output_json)


if __name__ == '__main__':
    cli()
