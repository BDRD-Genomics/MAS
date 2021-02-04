from django.core.management.base import BaseCommand
from django.db import connection


class Command(BaseCommand):
    help = "Database initialization SQL command to convert result_viewer_pdb_accession_mapping table to use " \
           "utf_bin collation (Allows for case-sensitivity in unique fields)"

    def handle(self, *args, **options):
        with connection.cursor() as cursor:
            cursor.execute("ALTER TABLE result_viewer_pdb_accession_mapping CONVERT TO CHARACTER SET utf8 COLLATE utf8_bin")
