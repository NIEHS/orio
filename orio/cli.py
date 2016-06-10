#!/usr/bin/env python

"""
Command-line interface for data validation:


python cli.py feature_list \
    ../../../public/media/unt1hr.obsTSS.bed \
    ./validators/data/hg19.chromSizes.txt

python cli.py sort_vector \
    ../../../public/media/unt1hr.obsTSS.bed \
    ../../../public/media/wgEncodeBroadHistoneA549CtcfEtoh02Sig.sortVector.txt

python cli.py bigwig \
    ../../../data/encode/hg19/wgEncodeBroadHistone/wgEncodeBroadHistoneK562Cbpsc369Sig.bigWig \
    ./validators/data/hg19.chromSizes.txt


To pass negative values, you must specify as below, specify as an option; you'd
write -- -100 for example [docs](http://click.pocoo.org/6/arguments/#option-like-arguments).

python cli.py analysis \
    start -- -2500 50 100 \
    ../../../public/media/unt1hr.obsTSS.bed \
    ./validators/data/hg19.chromSizes.txt
"""

import click
import validators


def write_errs(errs):
    if errs:
        click.echo(click.style(errs, bg='red', fg='white'), err=True)


@click.group()
def validate():
    pass


@validate.command()
@click.argument('feature_list')
@click.argument('chrom_sizes_file')
def feature_list(feature_list, chrom_sizes_file):
    validator = validators.FeatureListValidator(feature_list, chrom_sizes_file)
    validator.validate()
    write_errs(validator.display_errors())


@validate.command()
@click.argument('feature_bed')
@click.argument('sort_vector')
def sort_vector(feature_bed, sort_vector):
    validator = validators.SortVectorValidator(feature_bed, sort_vector)
    validator.validate()
    write_errs(validator.display_errors())


@validate.command()
@click.argument('bigwig')
@click.argument('chrom_sizes_file')
def bigwig(bigwig, chrom_sizes_file):
    validator = validators.BigWigValidator(bigwig, chrom_sizes_file)
    validator.validate()
    write_errs(validator.display_errors())


@validate.command()
@click.argument('bin_anchor',
                type=click.Choice(validators.AnalysisValidator.ANCHOR_OPTIONS))
@click.argument('bin_start', type=int)
@click.argument('bin_number', type=int)
@click.argument('bin_size', type=int)
@click.argument('feature_bed', type=str)
@click.argument('chrom_sizes', type=str)
@click.option('--stranded_bed', is_flag=True, help='Expect stranded bed')
def analysis(bin_anchor, bin_start, bin_number,
             bin_size, feature_bed, chrom_sizes,
             stranded_bed):
    validator = validators.AnalysisValidator(
        bin_anchor, bin_start, bin_number,
        bin_size, feature_bed, chrom_sizes,
        stranded_bed)
    validator.validate()
    write_errs(validator.display_errors())


if __name__ == '__main__':
    validate()
