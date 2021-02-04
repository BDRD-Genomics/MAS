#!/usr/bin/env python
'''
Matt Lueder

This script is for uploading phages to our internal database. To use this script, simply pass the complete fasta file
to it as an argument. The name which is entered in the database is taken from the fasta filename, so make sure this
is correct!

To see usage information use the option '-h'.
'''
import sys, os, django

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "LIMS.settings")
django.setup()

from genome import models as amd_models
import argparse
import re
import sys
import os
from Bio import SeqIO

from amd_database_scripts.paths import *

def parse_args():
    arg_parser = argparse.ArgumentParser(
        description='This script is for uploading phages to our internal database. To use this script, simply pass the '
                    'complete fasta file to it as an argument. The name which is entered in the database is taken from '
                    'the fasta filename, so make sure this is correct!')

    arg_parser.add_argument(
        '-f', dest='force', action='store_true',
        help='Overwrite previous entries and bypass name check.')

    arg_parser.add_argument('fasta',
                            help='Fasta file containing the complete phage genome to be uploaded. File should be named'
                                 ' the same as the phage.')

    return arg_parser.parse_args()


def test_phage_name(name):
    '''
    Test phage name to see if it is in the standard format. If not print warning and ask for user input_fs to continue.
    '''
    pattern = re.compile(r"AMD_[A-Z]_[a-z]+_[0-9A-Z]+_Phi_[0-9]{3}$")
    if pattern.match(name) is None:
        print("Genome name %s is not in the standard format. Upload anyways?" % name)
        while True:
            response = input('(Y/n) ').lower()
            if response in ['y', 'yes']:
                return

            elif response in ['n', 'no']:
                sys.exit(0)


def main(args):
    '''
    Start script execution
    :param args: Namespace object from argparse
    '''
    # Parse name of phage from filename
    file_name = os.path.split(args.fasta)[1]
    phage_name = file_name[:file_name.find('.')]

    if not args.force:
        test_phage_name(phage_name)

    # Read genome from file
    genome = SeqIO.read(args.fasta, 'fasta').seq.__str__().upper()

    # Connect to database
    # connection = sqlite3.connect(DB_FILE)
    # cursor = connection.cursor()

    # try:
    obj = amd_models.Phages(name=phage_name, genome=genome)
        # cursor.execute('INSERT INTO phages VALUES(?, ?);', (phage_name, genome))


    # except sqlite3.IntegrityError as e:
    #     # If protein already annotated, we update only if -f arg was supplied
    #     if "UNIQUE constraint failed: phages.name" in e.__str__() and args.force:
    #         cursor.execute('UPDATE phages SET genome = ? WHERE name = ?;', (genome, phage_name))
    #
    #     else:
    #         print(e)

    # connection.commit()
    # connection.close()
    obj.save()


if __name__ == '__main__':
    main(parse_args())

