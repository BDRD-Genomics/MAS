#!/usr/bin/env python
'''
Matthew Lueder

Creates BLAST database(s). Can create a create a nucleotide db containing all complete genomes from the phages table,
and a protein db with all annotated proteins from the annotations table.

To see usage information use the option '-h'.
'''
import sqlite3
import argparse
import os
import subprocess
from tempfile import TemporaryDirectory

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import create_deliverables


def parse_args():
    arg_parser = argparse.ArgumentParser(
        description='Creates BLAST database(s). Can create a create a nucleotide db containing all complete genomes '
                    'from the phages table, and a protein db with all polished proteins from the annotations table.')

    arg_parser.add_argument(
        '-n', dest='create_nucl_db', action='store_true',
        help='Create a nucleotide BLAST db with all complete genomes from the phages table')

    arg_parser.add_argument(
        '-p', dest='create_protein_db', action='store_true',
        help='Create a protein BLAST db with all polished proteins from the annotations table')

    arg_parser.add_argument('output_folder', type=str, help='Destination of output files.')

    return arg_parser.parse_args()


def main(args):
    '''
    Start script execution
    :param args: Namespace object from argparse
    '''
    # Connect to database
    connection = sqlite3.connect(DB_FILE)
    cursor = connection.cursor()

    result = cursor.execute('SELECT * FROM phages')
    phages = result.fetchall()

    if args.create_nucl_db:
        records = []

        for phage in phages:
            name = phage[0]
            # genome = Seq(phage[1])
            genome = Seq(phage[1], IUPAC.unambiguous_dna)

            record = SeqRecord(genome, id=name, description='')
            records.append(record)

        nucl_fa = os.path.join(args.output_folder, 'AMD_Complete_Phages.fa')
        SeqIO.write(records, nucl_fa, 'fasta')

        # Use the fasta file to create a blast db with makeblastdb
        nucl_db = os.path.join(args.output_folder, 'AMD_Complete_Phages')
        subprocess.run(['makeblastdb', '-in', '"%s"' % nucl_fa, '-input_type', 'fasta', '-dbtype', 'nucl',
                        '-title', 'MAS Genomes', '-out', nucl_db], check=True)

    if args.create_protein_db:
        # Save faa files for each individual phage to a temp directory
        with TemporaryDirectory() as temp_dir:
            paths = []

            for phage in phages:
                name = phage[0]

                cd_args = argparse.Namespace(
                    output_folder=temp_dir,
                    phage_name=name,
                    deliverables=['faa'],
                    only_annotated=True)

                try:
                    create_deliverables.main(cd_args)
                    paths.append(os.path.join(temp_dir, '%s.faa' % name))

                except Exception as e:
                    print(e)

            # Concat all faa files in temp dir to a fasta in the output dir
            annotated_proteins_path = os.path.join(args.output_folder, 'AMD_Annotated_Proteins.faa')
            subprocess.run('cat ' + os.path.join(temp_dir, '*'), shell=True, stdout=open(annotated_proteins_path, 'w+'))

            # Use the fasta file to create a blast db with makeblastdb
            prot_db = os.path.join(args.output_folder, 'AMD_Annotated_Proteins')
            subprocess.run(['makeblastdb', '-in', '"%s"' % annotated_proteins_path, '-input_type', 'fasta', '-dbtype',
                            'prot', '-title', 'MAS Proteins', '-out', prot_db], check=True)


if __name__ == '__main__':
    main(parse_args())