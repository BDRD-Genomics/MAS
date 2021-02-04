from django.db import models
from django.contrib.auth.models import User
from simple_history.models import HistoricalRecords

from genome.models import *

import os


def result_upload_to(instance, filename):
    return os.path.join('search_results', instance.__class__.__name__, instance.database, filename)


search_status_options = (
    (0, 'complete'),
    (1, 'running'),
    (2, 'error')
)


class UserPreferences(models.Model):
    user = models.OneToOneField(User, on_delete=models.CASCADE)
    dark_mode_activated = models.BooleanField(default=False)


class HHSearch_Result(models.Model):
    database_options = (
        ('pdb', 'Protein Data Bank'),
    )

    annotation = models.ForeignKey(Annotation, on_delete=models.CASCADE, null=False, db_index=True)
    database = models.CharField(max_length=50, choices=database_options, null=False, blank=False)
    result = models.FileField(upload_to=result_upload_to, null=True)
    run_date = models.DateTimeField(null=True)
    status = models.IntegerField(default=0, choices=search_status_options)

    class Meta:
        # db_table = 'genome_hhsearch_result'
        unique_together = ['annotation', 'database']


class Blastp_Result(models.Model):
    database_options = (
        ('swissprot', 'SwissProt'),
        ('internal', 'Previously Annotated Proteins'),
        ('nr', 'NCBI Protein'),
    )

    annotation = models.ForeignKey(Annotation, on_delete=models.CASCADE, null=False, db_index=True)
    database = models.CharField(max_length=50, choices=database_options, null=False, blank=False)
    result = models.FileField(upload_to=result_upload_to, null=True)
    run_date = models.DateTimeField(null=True)
    status = models.IntegerField(default=0, choices=search_status_options)

    class Meta:
        # db_table = 'genome_blastp_result'
        unique_together = ['annotation', 'database']


class RPSBlast_Result(models.Model):
    database_options = (
        ('cdd', 'Conserved Domain Database'),
    )

    annotation = models.ForeignKey(Annotation, on_delete=models.CASCADE, null=False, db_index=True)
    database = models.CharField(max_length=50, choices=database_options, null=False, blank=False)
    result = models.FileField(upload_to=result_upload_to, null=True)
    run_date = models.DateTimeField(null=True)
    status = models.IntegerField(default=0, choices=search_status_options)

    class Meta:
        # db_table = 'genome_rpsblast_result'
        unique_together = ['annotation', 'database']


class PDB_Accession_Mapping(models.Model):
    pdb_accession = models.CharField(max_length=20, null=False, blank=False, unique=True)
    pdb_chain_name = models.CharField(max_length=1500, blank=False, null=True)
