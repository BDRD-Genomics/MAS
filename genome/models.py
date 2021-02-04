from django.db import models
import re
from django.core.exceptions import ValidationError
from django.shortcuts import redirect
from django.contrib.auth.models import User
from simple_history.models import HistoricalRecords
from django.conf import settings
import django.dispatch


# Validators
def validate_genome_name(value):
    if settings.GENOME_NAME_FORMAT:
        pattern = re.compile(settings.GENOME_NAME_FORMAT)
        if pattern.match(value) is None:
            raise ValidationError("%s is not a valid phage name. Would You like to continue with this name?" % value, code=1)


def validate_phage_error(option):
    if option != 'Yes':
        raise ValidationError("Genome was not accepted", code=2)
    else:
        redirect('genome:phage_list')


def validate_duplicate_name(name):
    # phages = Features.objects.filter(name=name)
    genomes = Genome.objects.all()
    for genome in genomes:
        if genome.genome_name == name:
            raise ValidationError("%s is already a genome name." % name, code=3)


genome_upload_complete = django.dispatch.Signal()


# Create your models here.
class Genome(models.Model):
    genome_name = models.CharField(max_length=100, unique=True)
    genome_sequence = models.TextField(max_length=15000000)
    organism = models.CharField(max_length=100, default='phage')

    def __str__(self):
        return self.genome_name


class Feature(models.Model):
    feature_options = (
        ('CDS', 'Coding Sequence'),
        ('Repeat Region', 'Repeat Region'),
        ('tRNA', 'tRNA'),
    )

    """
    on_delete=models.CASCADE() - deletes features of a phage when a phage is deleted
    """
    genome = models.ForeignKey('Genome', on_delete=models.CASCADE)

    start = models.IntegerField(default=0)
    stop = models.IntegerField(default=0)

    #
    type = models.CharField(max_length=50, choices=feature_options)

    # strand = '+' or '-' depending on backbone sidegenomes_phage
    strand = models.CharField(max_length=1)

    # obtain information from the phage annotations
    annotation = models.ForeignKey('Annotation', blank=True, null=True, on_delete=models.PROTECT)

    def __str__(self):
        return "%s: %s %s..%s %s" % (self.type, self.genome, self.start, self.stop, self.strand)


class Annotation(models.Model):
    # if you update this, you need to update the flag dict in confirm_upload_annotations in genome views
    flag_options = (
        (0, 'GREEN'),
        (1, 'YELLOW'),
        (2, 'RED'),
        (3, 'REVIEW NAME'),
        (4, 'N/A'),
        (5, 'ORANGE'),
        (6, 'ENDOLYSIN'),
        (7, 'UNANNOTATED'),
        (8, 'tRNA'),
        (9, 'TERMINAL REPEAT')
    )
    annotation = models.CharField(max_length=255, blank=True, null=True, default='No Annotation')
    #Amino Acid Sequence
    sequence = models.TextField(max_length=10000, unique=True)
    public_notes = models.TextField(max_length=30000, blank=True, null=True, default='')
    private_notes = models.TextField(max_length=30000, blank=True, null=True, default='')
    # default 7 for unannotated phage
    flag = models.IntegerField(default=7, choices=flag_options)
    assigned_to = models.ForeignKey(User, on_delete=models.SET_NULL, null=True)
    history = HistoricalRecords()

    @property
    def accession(self):
        chars = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        integer = abs(self.id)
        result = ''
        padding = []

        while integer > 0:
            integer, remainder = divmod(integer, 36)
            padding.insert(0, chars[remainder])
        if len(padding) < 5:
            while len(padding) < 5:
                padding.insert(0, '0')
        for i in padding:
            result = result + i
        return result

    def __str__(self):
        return "%s | %s" % (self.annotation, self.get_flag_display())

