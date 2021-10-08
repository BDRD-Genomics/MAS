import subprocess
import sys
from tempfile import TemporaryDirectory
import os, django
from argparse import Namespace

from genome.genomic_loci_conversions import *

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


def parse_glimmer_results(path):
    with open(path, 'r') as glimmer_predict:
        for line in glimmer_predict:
            if line.startswith(">"):
                continue
            else:
                lineSplit = [x.strip() for x in line.split()]

                start, stop, strand = glimmer_to_db_standard(
                    start=int(lineSplit[1]),
                    stop=int(lineSplit[2]),
                    strand=lineSplit[3][0]
                )

                yield Namespace(
                    start=start,
                    stop=stop,
                    strand=strand
                )


def parse_trnascan_results(path):
    with open(path, 'r') as trnascan_results:
        ready = False
        for line in trnascan_results:
            if line.startswith('---'):
                ready = True
                continue

            if ready is True:
                lineSplit = [x.strip() for x in line.split()]

                begin = int(lineSplit[2])
                end = int(lineSplit[3])

                # Convert to database standard genomic loci representation
                start, stop, strand = trnascan_to_db_standard(begin, end)

                yield Namespace(
                    start=start,
                    stop=stop,
                    strand=strand,
                    amino=lineSplit[4],
                    codon=lineSplit[5],
                    score=float(lineSplit[8])
                )