#!/usr/bin/env python
'''
Matt Lueder

This script uploads a phage genome, and uses the genome as input into glimmer3 and trnascan-se to automatically add
cds and trna features. It wraps the upload_phage.py and create_features_for_phage.py scripts. Annotations for tRNAs
and terminal repeats are created automatically.

Requires the tools long-orfs, extract, build-icm, glimmer3, and trnascan-se to be in PATH.
These are available as conda packages:
> conda install glimmer trnascan-se

To see usage information use the option '-h'.
'''

import subprocess
import sys
from tempfile import TemporaryDirectory
import os, django

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "LIMS.settings")
django.setup()


def run_glimmer(fasta_path, phage_name, output_destination, minimum_gene_len=90):
    '''
    Run glimmer and return path to predict file.
    '''
    print('*** RUNNING GLIMMER ***')

    # Run long-orfs
    long_orf_file = open(os.path.join(output_destination, '%s.glimmer.longorfs' % phage_name), 'w')
    subprocess.run(['long-orfs', '-n', '-t', '1.1', fasta_path, '-'], stdout=long_orf_file, check=True)

    # Run extract
    extract_file = open(os.path.join(output_destination, '%s.glimmer.train' % phage_name), 'w+')
    subprocess.run(['extract', '-t', fasta_path, long_orf_file.name], stdout=extract_file, check=True)
    extract_file.seek(0)

    # Run build icm
    icm_file_path = os.path.join(output_destination, '%s.glimmer.icm' % phage_name)
    subprocess.run(['build-icm', '-r', icm_file_path], check=True, stdin=extract_file)

    # Run glimmer3
    subprocess.run(['glimmer3', '-o50', '-g{}'.format(minimum_gene_len), '-t30', fasta_path, icm_file_path,
                    os.path.join(output_destination, '%s.glimmer' % phage_name)], check=True)

    long_orf_file.close()
    extract_file.close()

    return os.path.join(output_destination, '%s.glimmer.predict' % phage_name)


def run_trnascan_se(fasta_path, phage_name, output_destination):
    '''
    Run tRNAscan-SE and return output file
    '''
    print('*** RUNNING GLIMMER ***')

    out_file_path = os.path.join(output_destination, '%s.trnascan_se.results.txt' % phage_name)
    subprocess.run(['tRNAscan-SE', '-B', '-D', '-o', out_file_path, fasta_path], check=True)
    return out_file_path
