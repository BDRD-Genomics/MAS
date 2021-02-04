import pandas as pd

from amd_database_scripts.genomic_loci_conversions import *
from django.core.exceptions import ObjectDoesNotExist
from Bio.Blast import NCBIXML

from result_viewer.models import Feature, Blastp_Result


def blastp_alignment_to_dict(alignment):
    '''
    :param hsp: biopython Alignment object from NCBIXML parsing (Don't use SearchIO)
    :return: A dict which can be used in jinja templating and serialized to json
    '''
    hits = []
    top_score = 0
    sum_score = 0
    lowest_evalue = 100
    query_start = None
    query_end = None

    for hsp in alignment.hsps:
        hits.append({
            'hit_start': hsp.sbjct_start,
            'hit_end': hsp.sbjct_end,
            'query_start': hsp.query_start,
            'query_end': hsp.query_end,
            'evalue': hsp.expect,
            'score': hsp.bits
        })
        lowest_evalue = hsp.expect if hsp.expect < lowest_evalue else lowest_evalue
        top_score = hsp.bits if hsp.bits > top_score else top_score
        query_start = hsp.query_start if query_start is None or hsp.query_start < query_start else query_start
        query_end = hsp.query_end if query_end is None or hsp.query_end > query_end else query_end
        sum_score += hsp.bits

    return {
        'hit_id': alignment.accession,
        'hit_description': alignment.title,
        'hit_seq_len': alignment.length,
        'top_score': top_score,
        'sum_score': sum_score,
        'lowest_evalue': lowest_evalue,
        'query_start': query_start,
        'query_end': query_end,
        'hits': hits,
    }


def GetFeatureTable(phage_name, genome):
    '''
    Return a Pandas DataFrame which contains information about a phage's features.
    :param phage_name: (str) Name of phage
    :param genome: (Bio.Seq) Full genome sequence
    :return: Pandas DataFrame
    '''
    # Connect to database
    # connection = sqlite3.connect(DB_FILE)
    # cursor = connection.cursor()
    #
    # # Get phage's features
    # result = cursor.execute('SELECT * FROM features WHERE phage_name = ?', (phage_name,))
    # features = result.fetchall()
    # if len(features) == 0:
    #     raise RuntimeError('ERROR: %s does not have any features.' % phage_name)

    # Get all features associated with phage
    features = []
    feature_obj_list = list(Feature.objects.filter(genome__genome_name=phage_name))
    for feature_obj in feature_obj_list:
        annotation_obj = feature_obj.annotation

        dna = get_dna_sequence(feature_obj.start, feature_obj.stop, feature_obj.strand, genome)
        blastp_result_string = ''
        if feature_obj.type == 'CDS':
            prot = get_protein_sequence(feature_obj.start, feature_obj.stop, feature_obj.strand, genome)
            # try:
            #     #TODO: For downloading excel files for telework
            #     db = 'internal'
            #     blastp_result = Blastp_Result.objects.get(annotation=annotation_obj, database=db)
            #     if blastp_result.result:
            #         record = NCBIXML.read(blastp_result.result.open(mode='r'))
            #         alignments = [blastp_alignment_to_dict(x) for x in record.alignments]
            #
            #         for alignment in alignments:
            #             if blastp_result_string:
            #                 blastp_result_string = '%s; %s' % (blastp_result_string, alignment['hit_description'])
            #             else:
            #                 blastp_result_string = alignment['hit_description']
            #
            # except ObjectDoesNotExist:
            #     pass
        else:
            prot = ''

        features.append((
            phage_name,
            feature_obj.start,
            feature_obj.stop,
            feature_obj.strand,
            feature_obj.type,
            len(dna),
            len(prot),
            annotation_obj.annotation,
            annotation_obj.public_notes,
            annotation_obj.private_notes,
            annotation_obj.get_flag_display(),
            dna.__str__(),
            prot.__str__(),
            blastp_result_string,
        ))

    # Store data in a dataframe
    df = pd.DataFrame(data=features,
                      columns=['Genome Name', 'Feature Start', 'Feature End', 'Strand', 'Feature Type', 'Gene Length',
                               'Protein Length', 'Manual Annotation', 'Public Notes', 'Internal Notes', 'Flag',
                               'Gene Sequence', 'Protein Sequence', 'Internal BLASTp Hits']
                      )

    # Sort by feature starting location
    df = df.sort_values('Feature Start', ascending=True)

    return df
