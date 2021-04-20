from result_viewer.models import *
from genome.models import *

from django.core.exceptions import ObjectDoesNotExist
from django.contrib.auth.models import User

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from simple_history.utils import get_history_manager_for_model


def move_users():
    print('Starting moving users.')
    for user in User.objects.using('old-mas').iterator():
        user.id = None
        if len(User.objects.using('default').filter(username=user.username)) == 0:
            user.save(using='default')
    print('Done moving users.')


def move_annotations():
    '''
    1) Move all annotations
    2) Move phage one-by-one
    3) Move features associated with phage
    4) Move search results
    '''
    print('Starting moving annotations.')
    histories = []
    annos = []
    useqs = set()
    for annotation in Annotation.objects.using('old-mas').all():
        if annotation.sequence in useqs:
            continue
        else:
            useqs.add(annotation.sequence)

        for history in annotation.history.using('old-mas').iterator():
            history.history_id = None

            if history.assigned_to_id:
                assigned_to = User.objects.using('old-mas').get(id=history.assigned_to_id)
                history.assigned_to = User.objects.get(username=assigned_to.username)

            if history.history_user_id:
                history_user = User.objects.using('old-mas').get(id=history.history_user_id)
                history.history_user = User.objects.get(username=history_user.username)
            # history.save(using='default')
            histories.append(history)

        annotation.id = None
        if annotation.assigned_to_id:
            assigned_to = User.objects.using('old-mas').get(id=annotation.assigned_to_id)
            annotation.assigned_to = User.objects.get(username=assigned_to.username)
        annos.append(annotation)
        # annotation.save(using='default')

    Annotation.objects.using('default').bulk_create(annos, batch_size=1000)
    get_history_manager_for_model(Annotation).using('default').bulk_create(histories, batch_size=1000)
    print('Done moving Annotations.')


def move_phage_and_features():
    print('Starting moving phages and features.')
    features = []
    new_annotations = []

    i = 1
    for phage in Phage.objects.using('old-mas').all():
        print('Starting phage {}'.format(i))
        i += 1

        genome = Genome(
            genome_name=phage.phage_name,
            genome_sequence=phage.genome
        )
        genome.save(using='default')

        phage_features = old_Feature.objects.using('old-mas').filter(phage__phage_name=phage.phage_name)

        for feature in phage_features.all():
            # Find annotation associated with feature
            ftype = None
            if feature.type == 'CDS':
                sequence = get_protein_sequence(
                    feature.start,
                    feature.stop,
                    feature.strand,
                    Seq(genome.genome_sequence, IUPAC.unambiguous_dna)
                )
                ftype = 'CDS'

            elif feature.type == 'repeat_region' or feature.type == 'Repeat Region':
                sequence = get_dna_sequence(
                    feature.start,
                    feature.stop,
                    feature.strand,
                    Seq(genome.genome_sequence, IUPAC.unambiguous_dna)
                )
                ftype = 'Repeat Region'

            elif feature.type == 'tRNA':
                sequence = get_rna_sequence(
                    feature.start,
                    feature.stop,
                    feature.strand,
                    Seq(genome.genome_sequence, IUPAC.unambiguous_dna)
                )
                ftype = 'tRNA'

            else:
                raise ValueError('Feature type not recognized')

            try:
                # If not found in amd table, look if we already created it for MAS table
                mas_annotation = Annotation.objects.using('default').get(sequence=sequence)

            # start from here
            except ObjectDoesNotExist:
                print('Warning: Annotation not found. Creating...')
                mas_annotation = Annotation(
                    sequence=sequence,
                    flag=7
                )
                new_annotations.append(mas_annotation)
                mas_annotation.save()

            mas_feature = Feature(
                genome=genome,
                start=feature.start,
                stop=feature.stop,
                type=ftype,
                strand=feature.strand,
                annotation=mas_annotation
            )
            features.append(mas_feature)

    Feature.objects.using('default').bulk_create(features, batch_size=1000)
    print('Done moving phages and features')


def move_search_results():
    print('Starting moving search results.')

    rpsblast_bin = []
    hhsearch_bin = []
    blastp_bin = []
    for old_search_model in [old_Blastp_Result, old_HHSearch_Result, old_RPSBlast_Result]:
        for old_search_result in old_search_model.objects.using('old-mas').all():
            # Get annotation this search result links to
            old_anno = Annotation.objects.using('old-mas').get(id=old_search_result.annotation_id)
            mas_annotation = Annotation.objects.using('default').get(sequence=old_anno.sequence)

            if isinstance(old_search_result, old_HHSearch_Result):
                result_class = HHSearch_Result
                bin = hhsearch_bin
            elif isinstance(old_search_result, old_Blastp_Result):
                result_class = Blastp_Result
                bin = blastp_bin
            elif isinstance(old_search_result, old_RPSBlast_Result):
                result_class = RPSBlast_Result
                bin = rpsblast_bin
            else:
                raise RuntimeError('Invalid class: "{}"'.format(old_search_result.__class__))

            obj = result_class(
                annotation=mas_annotation,
                database=old_search_result.database,
                run_date=old_search_result.run_date,
                status=old_search_result.status,
                result=old_search_result.result
            )
            # results.append(obj)
            # obj.result.save(name=base36encode(mas_annotation.id)
            # obj.save(using='default')
            bin.append(obj)

            # with io.StringIO(amd_search_result.hhpred_pdb) as fh:
            #     obj = HHSearch_Result(
            #         annotation=mas_annotation,
            #         database='pdb',
            #         run_date=datetime(2019, 5, 9)
            #     )
            #     obj.result.save(name=base36encode(mas_annotation.id), content=ContentFile(content=fh.read()))
            #     obj.save()
        print('Done binning {}'.format(old_search_model))
    RPSBlast_Result.objects.using('default').bulk_create(rpsblast_bin, batch_size=1000)
    HHSearch_Result.objects.using('default').bulk_create(hhsearch_bin, batch_size=1000)
    Blastp_Result.objects.using('default').bulk_create(blastp_bin, batch_size=1000)
    print('Done moving search results.')


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