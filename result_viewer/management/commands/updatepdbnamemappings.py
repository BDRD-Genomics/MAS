from django.core.management.base import BaseCommand, CommandError
import sys

from result_viewer.models import *

from Bio import SeqIO


class Command(BaseCommand):
    help = 'Fills database with PDB accession -> PDB chain name mapping. Will add mapping to uniprot data later.'

    def add_arguments(self, parser):
        parser.add_argument('pdb_seqres_fasta', type=str,
                            help='Can be found here: ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt')
        parser.add_argument('--override_existing', action='store_true',
                            help='If mapping data exists in database, delete it and use provided file to create a new '
                                 'mapping. Otherwise, if mapping data exists, do nothing.')

    def handle(self, *args, **options):
        with open(options['pdb_seqres_fasta'], 'r') as pdb_seqres_fasta:
            if PDB_Accession_Mapping.objects.all().count() > 0 and not options['override_existing']:
                sys.exit(0)

            # Clear all old data
            PDB_Accession_Mapping.objects.all().delete()

            # Load in PDB -> Chain name mappings
            new_objects = []
            for record in SeqIO.parse(pdb_seqres_fasta, 'fasta'):
                accession = record.id[:4].upper() + record.id[4:]
                name = record.description.split('  ')[-1]

                obj = PDB_Accession_Mapping(pdb_accession=accession, pdb_chain_name=name)
                new_objects.append(obj)

            PDB_Accession_Mapping.objects.bulk_create(new_objects, batch_size=2500)
