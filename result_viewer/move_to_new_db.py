from result_viewer.models import *

from django.core.exceptions import ObjectDoesNotExist
from django.core.files.base import ContentFile
from django.db.utils import IntegrityError

import io
from datetime import datetime

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def move_annotations():
    '''
    1) Move all annotations
    2) Move phage one-by-one
    3) Move features associated with phage
    4) Move search results
    '''
    for amd_annotation in Annotations.objects.iterator():

        if amd_annotation.annotation.startswith('tRNA-') and amd_annotation.flag == 4:
            mas_flag = 8
        elif amd_annotation.gbk_notes.startswith('Direct terminal repeat.') and amd_annotation.flag == 4:
            mas_flag = 9
        else:
            mas_flag = amd_annotation.flag

        try:
            mas_annotation = Annotation(
                sequence=amd_annotation.sequence,
                annotation=amd_annotation.annotation,
                public_notes=amd_annotation.gbk_notes,
                private_notes=amd_annotation.private_notes,
                flag=mas_flag
            )
            mas_annotation.save()

        except IntegrityError:
            Annotation.objects.filter(sequence=amd_annotation.sequence).update(
                annotation=amd_annotation.annotation,
                public_notes=amd_annotation.gbk_notes,
                private_notes=amd_annotation.private_notes,
                flag=mas_flag
            )


phages = []
phage_features = []  # each element a list of features which belong to phage at same index
new_annotations = []

def move_phage_and_features():
    phages = []
    phage_features = []  # each element a list of features which belong to phage at same index
    new_annotations = []

    try:
        for amd_phage in Phages.objects.iterator():
            mas_phage = Genome(
                phage_name=amd_phage.name,
                genome=amd_phage.genome
            )
            phages.append(mas_phage)
            mas_phage.save()

            amd_features = Features.objects.filter(phage_name=amd_phage.name)

            features = []  # The current phage's features (to be appended to phage_features list)
            for amd_feature in amd_features.iterator():
                # Find annotation associated with feature
                if amd_feature.type == 'CDS':
                    sequence = get_protein_sequence(
                        amd_feature.start,
                        amd_feature.end,
                        amd_feature.strand,
                        Seq(amd_phage.genome, IUPAC.unambiguous_dna)
                    )

                elif amd_feature.type == 'repeat_region':
                    sequence = get_dna_sequence(
                        amd_feature.start,
                        amd_feature.end,
                        amd_feature.strand,
                        Seq(amd_phage.genome, IUPAC.unambiguous_dna)
                    )

                elif amd_feature.type == 'tRNA':
                    sequence = get_rna_sequence(
                        amd_feature.start,
                        amd_feature.end,
                        amd_feature.strand,
                        Seq(amd_phage.genome, IUPAC.unambiguous_dna)
                    )

                else:
                    raise ValueError('Feature type not recognized')

                try:
                    # If not found in amd table, look if we already created it for MAS table
                    mas_annotation = Annotation.objects.get(sequence=sequence)

                except ObjectDoesNotExist:
                    mas_annotation = Annotation(
                        sequence=sequence,
                        flag=7
                    )
                    new_annotations.append(mas_annotation)
                    mas_annotation.save()

                mas_feature = Feature(
                    phage=mas_phage,
                    start=amd_feature.start,
                    stop=amd_feature.end,
                    type=amd_feature.type,
                    strand=amd_feature.strand,
                    annotation=mas_annotation
                )
                features.append(mas_feature)

            phage_features.append(features)

        for i, phage in enumerate(phages):
            #     print(phage.name)
            #     phage.save()
            for feature in phage_features[i]:
                print(feature)
                feature.save()

    except:
        for annoation in new_annotations:
            annoation.delete()

        for i, phage in enumerate(phages):
            phage.delete()

            for feature in phage_features[i]:
                feature.delete()


def move_search_results():
    for amd_search_result in SearchResults.objects.iterator():
        # Get annotation this search result links to
        mas_annotation = Annotation.objects.get(sequence=amd_search_result.sequence)

        if amd_search_result.hhpred_pdb:
            with io.StringIO(amd_search_result.hhpred_pdb) as fh:
                obj = HHSearch_Result(
                    annotation=mas_annotation,
                    database='pdb',
                    run_date=datetime(2019, 5, 9)
                )
                obj.result.save(name=base36encode(mas_annotation.id), content=ContentFile(content=fh.read()))
                obj.save()

        if amd_search_result.blastp_swissprot:
            with io.StringIO(amd_search_result.blastp_swissprot) as fh:
                obj = Blastp_Result(
                    annotation=mas_annotation,
                    database='swissprot',
                    run_date=datetime(2019, 5, 9)
                )
                obj.result.save(name=base36encode(mas_annotation.id), content=ContentFile(content=fh.read()))
                obj.save()

        if amd_search_result.rpsblast_cdd:
            with io.StringIO(amd_search_result.rpsblast_cdd) as fh:
                obj = RPSBlast_Result(
                    annotation=mas_annotation,
                    database='cdd',
                    run_date=datetime(2019, 5, 9)
                )
                obj.result.save(name=base36encode(mas_annotation.id), content=ContentFile(content=fh.read()))
                obj.save()




'''
BEGINNING CODE FROM AMD DATABASE SCRIPTS - genomic_loci_conversions.py

Matthew Lueder

This is for converting descriptions of feature locations to match the indexing/format of various file formats.

The database standard format is:
  - 0-based (count starts from zero)
  - half-open (last base not included in feature)
  - 'strand' represented by '+' or '-'
  - start can > stop when feature is on negative or when feature stretches over termini
'''


def trnascan_to_db_standard(begin, end):
    '''
    tRNAscan-SE results are:
      - 1-based (count starts from one)
      - Closed (last base is included)
      - Can not represent tRNAs which stretch over ends of circular genome
      - Do not give a strand. Negative strand tRNAs have a start > stop
    '''
    strand = '+' if begin < end else '-'

    if strand == '+':
        begin -= 1

    elif strand == '-':
        end -= 1

    return begin, end, strand


def glimmer_to_db_standard(start, stop, strand):
    '''
    Glimmer results are:
      - 1-based
      - Closed (last base is included)
      - 'strand' represented by '+' or '-'
      - start can > stop when feature is on negative or when feature stretches over termini
      - includes stop codon
    '''
    if strand == '+':
        start -= 1

    elif strand == '-':
        stop -= 1

    else:
        raise ValueError("'strand' must be either '+' or '-'")

    return start, stop, strand


def db_standard_to_gff(start, stop, strand, genome_length):
    '''
    GFF version 3 results are:
      - 1-based
      - Closed (last base is included)
      - start site always greater than end (End greater than genome length for wrap-arounds)
    '''
    begin = min(start, stop) + 1
    end = max(start, stop)

    # Handle forward wrap arounds
    if strand == '+' and start > stop:
        begin = start + 1
        end = genome_length + stop

    return begin, end, strand


def db_standard_to_tbl_text(start, stop, strand, genome_length, feature_type):
    '''
    Outputs text version of coordinates for tbl files.
    Tbl results are:
      - 1-based
      - Closed (last base is included)
      - No strand - start > stop for features on - strand
      - wrap-around genes encoded with four coords
    '''
    if strand == '+':
        start += 1

    elif strand == '-':
        stop += 1

    if start < stop and strand == '-':
        return '%i\t1\t%s\n%i\t%i\n' % (start, feature_type, genome_length, stop)

    elif start > stop and strand == '+':
        return '%i\t%i\t%s\n1\t%i\n' % (start, genome_length, feature_type, stop)

    else:
        return '%i\t%i\t%s\n' % (start, stop, feature_type)


def get_dna_sequence(start, stop, strand, genome):
    if start < stop and strand == '+':
        return genome[start:stop]

    elif start < stop and strand == '-':
        return (genome[stop:] + genome[:start]).reverse_complement()

    elif start > stop and strand == '+':
        return (genome[start:] + genome[:stop])

    elif start > stop and strand == '-':
        return genome[stop:start].reverse_complement()


def get_protein_sequence(start, stop, strand, genome):
    '''
    Extract protein sequence from genome given coordinates in db standard format.
    'genome' is a Bio.Seq object
    '''
    return get_dna_sequence(start, stop, strand, genome).translate(table=11, cds=True)


def get_rna_sequence(start, stop, strand, genome):
    return get_dna_sequence(start, stop, strand, genome).transcribe()

# with padding
# 13 in base 10 = 000D in base 36
def base36encode(integer):
    chars = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer = abs(integer)
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


