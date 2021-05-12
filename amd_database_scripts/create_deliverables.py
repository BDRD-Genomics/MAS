#!/usr/bin/env python
'''
Matt Lueder

This script creates a fasta, genbank, gff, faa, tbl, sqn, and annotation table for the selected phage.

To see usage information use the option '-h'.
'''
# import sqlite3
import argparse
import subprocess
import sys
import re
import os
from datetime import datetime
import fileinput

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
import pandas as pd

from amd_database_scripts.genomic_loci_conversions import *
from amd_database_scripts.feature_table import GetFeatureTable

from result_viewer.models import Genome
from MAS import settings


def parse_args():
    arg_parser = argparse.ArgumentParser(
        description='This script creates a fasta, genbank, gff, faa, tbl, '
                    'sqn, and annotation table for the selected phage.')

    # arg_parser.add_argument('--only-annotated', dest='only_annotated', action='store_true',
    #                         help='Only include proteins which have an annotation in gbf, gff, faa, tbl, and sqn.')

    arg_parser.add_argument('output_folder', type=str, help='Destination of output files.')

    arg_parser.add_argument('phage_name', type=str, help='Name of the phage we would like to create deliverables for')

    arg_parser.add_argument('deliverables', type=str, nargs='*',
                            help='A list of deliverables to generate. If left empty, all deliverables will be created. '
                                 'Possible options are: fsa, faa, gff, tbl, gbf, xlsx, and sqn. '
                                 'Note: Generation of sqn requires the generation of the tbl and fsa. Generation of '
                                 'gbf requires generation of sqn, fsa, and tbl.')

    return arg_parser.parse_args()


def test_save(file):
    '''
    If file exists, ask user if they would like to overwrite it.
    :param file: Full file path to object we want to save
    :return: True if we should continue with save, false if not
    '''
    if os.path.exists(file):
        print("File with name of %s already exists in output path. Would you like to overwrite?" % file)
        while True:
            response = input('(Y/n) ').lower()
            if response in ['y', 'yes']:
                return True

            elif response in ['n', 'no']:
                return False

    else:
        return True


def write_fsa(genome, out_file, phage_name):
    '''
    Create a fsa (fasta) file representing the complete phage genome.
    Note: The 'fsa' extension is required by the tbl2asn tool
    :param genome: Biopython Seq object containing full phage genome
    :param out_file: file path to write to
    :param phage_name: Name of phage
    '''
    record = SeqRecord(genome, id=phage_name, description='')
    SeqIO.write(record, out_file, 'fasta')


def write_faa(df, out_file, phage_name, only_annotated=False):
    '''
    Create a faa file where each CDS is represented with a header and sequence
    :param df: Dataframe with all annotation data
    :param out_file: file path to write to
    :param phage_name: Name of phage
    :param only_annotated: Only output proteins which have an annotation
    '''
    groupings = df.groupby('Feature Type')['Genome Name'].count()
    # Skip writing of faa if there are no CDS
    if 'CDS' not in groupings:
        return
    num_cds = groupings['CDS']
    proteins = []
    prot_num = 0
    for i, row in df.iterrows():
        if row['Feature Type'] == 'CDS':
            prot_num += 1
            if only_annotated and row['Manual Annotation'] == '':
                continue

            protein_seq = Seq(row['Protein Sequence'], IUPAC.IUPACProtein)
            prot_id = '%s_%s' % (phage_name, str(prot_num).rjust(len(str(num_cds)), '0'))
            # Create string representation of locus
            loc = 'start=%i;end=%i;strand=%s' % (row['Feature Start'], row['Feature End'], row['Strand'])
            annotation = row['Manual Annotation'] if row['Manual Annotation'] != '' else 'Not yet annotated'
            pub_notes = row['Public Notes'] if row['Public Notes'] != '' else 'No public notes'
            internal_notes = row['Internal Notes'] if row['Internal Notes'] != '' else 'No internal notes'
            flag = row['Flag'] if row['Flag'] != '' else 'No Flag'
            description = '| %s | %s | %s | %s | %s' % (annotation, loc, pub_notes, internal_notes, flag)
            record = SeqRecord(protein_seq, id=prot_id, description=description)
            proteins.append(record)

    SeqIO.write(proteins, out_file, 'fasta')


def format_string_for_gff(string):
    '''
    Converts special characters in string so that the string can be used in a gff file
    '''
    conversions = {
        '\t': '%09',
        '\n': '%0A',
        '\r': '%0D',
        '%': '%25',
        ';': '%3B',
        '=': '%3D',
        '&': '%26',
        ',': '%2C'
    }
    conversions = {re.escape(k): conversions[k] for k in conversions}
    p = re.compile('|'.join(conversions.keys()))
    return p.sub(lambda m: conversions[re.escape(m.group(0))], string)


def write_gff(df, out_file, phage_name, genome_length):
    '''
    Create a gff which represents all features
    :param df: Dataframe with all annotation data
    :param out_file: file path to write to
    :param phage_name: Name of phage
    :param genome_length: Length of entire phage genome
    '''
    with open(out_file, 'w') as fh:
        # Print gff header
        fh.write('##gff-version 3\n')

        for i, row in df.iterrows():
            # Convert genomic coordinates to gff format
            start, stop, strand = db_standard_to_gff(row['Feature Start'], row['Feature End'],
                                                     row['Strand'], genome_length)

            # Get attributes column data
            if row['Feature Type'] == 'repeat_region' or row['Feature Type'] == 'Repeat Region':
                attributes = 'note=%s' % (format_string_for_gff(row['Public Notes']))

            elif row['Flag'] == 'UNANNOTATED':
                attributes = 'product=%s;' % format_string_for_gff('No Annotation')
                attributes = attributes + 'note=%s' % format_string_for_gff('')
            else:
                if row['Manual Annotation'] is None:
                    attributes = 'product=%s;' % format_string_for_gff('No Annotation')
                else:
                    attributes = 'product=%s;' % (format_string_for_gff(row['Manual Annotation']))
                if row['Public Notes'] is None:
                    attributes = attributes + 'note=%s' % format_string_for_gff('')
                else:
                    attributes = 'note=%s' % (format_string_for_gff(row['Public Notes']))

            # Write gff line
            line = '%s\t.\t%s\t%s\t%s\t.\t%s\t.\t%s\n' % (phage_name, row['Feature Type'],
                                                          start, stop, strand, attributes)
            fh.write(line)


def write_tbl(df, out_file, phage_name, genome_length):
    """
    Create a tbl file representing all features. This will be used to make the genbank file.
    :param df: Dataframe with all annotation data
    :param out_file: file path to write to
    :param phage_name: Name of phage
    :param genome_length: Length of entire phage genome
    :return:
    """
    with open(out_file, 'w') as fh:
        # Print tbl header
        # Print gff header
        fh.write('>Feature %s\n' % phage_name)

        for i, row in df.iterrows():
            # Write header of feature which describes location and feature type
            feature_header = db_standard_to_tbl_text(row['Feature Start'], row['Feature End'],
                                                     row['Strand'], genome_length, row['Feature Type'])
            fh.write(feature_header)

            # Write annotation
            if row['Feature Type'] != 'repeat_region':
                annotation = row['Manual Annotation'] if row['Manual Annotation'] != '' else 'Not yet annotated'
                fh.write('\t\t\tproduct\t%s\n' % annotation)

            # Write public notes
            if row['Public Notes'] is None:
                notes = ''
            else:
                notes = row['Public Notes'].replace('\n', '~')

            if notes != '':
                fh.write('\t\t\tnote\t%s\n' % notes)

            # Write trnaslation table
            fh.write('\t\t\ttransl_table\t11\n')


def collect_merge_values_for_docx(df, phage_name, genome):
    merge_values = {}

    merge_values['phage_name'] = phage_name
    merge_values['genome_size'] = '{:,}'.format(len(genome))
    merge_values['gc_content'] = '{:.2%}'.format(GC(genome) / 100)

    try:
        merge_values['host_name'] = phage_name[4:phage_name.find('_Phi_')]
    except:
        print('Warning: Could not fill host name field.', file=sys.stderr)

    num_cds, num_trna, num_hypothetical, dtr_length = [0] * 4
    notable_genes = ''

    for i, row in df.iterrows():
        if row['Feature Type'] == 'CDS':
            num_cds += 1
            if row['Manual Annotation'] is None:
                row['Manual Annotation'] = 'No Annotation'
            if 'hypothetical' in row['Manual Annotation'].lower():
                num_hypothetical += 1

            # output anything not flagged as GREEN, N/A, tRNA, or TERMINAL REPEAT to notable genes
            if row['Flag']:
                if row['Flag'] not in ['GREEN', 'N/A', 'tRNA', 'TERMINAL REPEAT', 'ENDOLYSIN']:
                    gene_info = '{}\n  - Flagged as {}'.format(row['Manual Annotation'], row['Flag'])
                    if row['Public Notes']:
                        gene_info = '{}\n  - {}'.format(gene_info, row['Public Notes'])
                    loc = 'Location: start=%i, end=%i, strand=%s' % \
                          (row['Feature Start'], row['Feature End'], row['Strand'])
                    gene_info = '{}\n  - {}'.format(gene_info, loc)

                    if notable_genes:
                        notable_genes = '{}\n\n{}'.format(notable_genes, gene_info)
                    else:
                        notable_genes = gene_info

        if row['Feature Type'] == 'tRNA':
            num_trna += 1

        if row['Feature Type'] == 'repeat_region':
            dtr_length = row['Feature End'] - row['Feature Start']

    merge_values['notable_genes'] = notable_genes if notable_genes else 'NSTR'
    merge_values['num_cds'] = '{:,}'.format(num_cds)
    merge_values['num_trna'] = '{:,}'.format(num_trna)
    merge_values['num_hypothetical'] = '{:,} ({:.2%} of proteins)'.format(num_hypothetical, num_hypothetical / num_cds)
    if dtr_length > 0:
        merge_values['dtr_length'] = '{:,}'.format(dtr_length)
    else:
        merge_values['dtr_length'] = 'N/A'

    merge_values['date_polished'] = datetime.today().strftime('%B %d, %Y')

    return merge_values


def write_docx(df, out_file, phage_name, genome):
    """
    Create a partially completed polished report for the phage. Some fields still need to be filled out manually.
    """
    from mailmerge import MailMerge

    merge_values = collect_merge_values_for_docx(df, phage_name, genome)

    # Create document
    template_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'polished_phage_report_template.docx')
    document = MailMerge(template_path)
    document.merge(**merge_values)
    document.write(out_file)


def shorten_phage_name(phage_name):
    host_gs = phage_name[4] + phage_name[6]
    phi_num = phage_name[phage_name.find('Phi_') + 4:phage_name.find('Phi_') + 7]
    shortened = host_gs + phi_num

    if phage_name.count('_', phage_name.find('Phi_')) > 1:
        shortened += phage_name[-1]

    return shortened

def main(args):
    '''
    Start script execution
    :param args: Namespace object from argparse
    '''
    # Make sure output destination exists
    if not os.path.isdir(args.output_folder):
        raise NotADirectoryError("Output directory '%s' does not exist." % args.output_folder)

    # Get phage's genome
    phage_obj = Genome.objects.get(genome_name=args.phage_name)
    genome = Seq(phage_obj.genome_sequence, IUPAC.ambiguous_dna)

    df = GetFeatureTable(args.phage_name, genome)

    if 'xlsx' in args.deliverables or len(args.deliverables) == 0:
        # Save the dataframe to an excel file
        excel_path = os.path.join(args.output_folder, '%s.xlsx' % args.phage_name)
        if test_save(excel_path):
            writer = pd.ExcelWriter(excel_path)
            df.to_excel(writer, 'Annotations')
            writer.save()

    if 'fsa' in args.deliverables or len(args.deliverables) == 0:
        # Write phage genome to fsa
        fsa_path = os.path.join(args.output_folder, '%s.fsa' % args.phage_name)
        if test_save(fsa_path):
            write_fsa(genome, fsa_path, args.phage_name)

    if 'faa' in args.deliverables or len(args.deliverables) == 0:
        # Save protein sequences to faa
        faa_path = os.path.join(args.output_folder, '%s.faa' % args.phage_name)
        if test_save(faa_path):
            # args.only_annotated is a hidden parameter which can be manually added to a Namespace object when calling
            # this script. It is used in create_blast_database.py
            if hasattr(args, 'only_annotated') and args.only_annotated == True:
                write_faa(df, faa_path, args.phage_name, args.only_annotated)
            else:
                write_faa(df, faa_path, args.phage_name)

    if 'gff' in args.deliverables or len(args.deliverables) == 0:
        # Write features to gff
        gff_path = os.path.join(args.output_folder, '%s.gff' % args.phage_name)
        if test_save(gff_path):
            write_gff(df, gff_path, args.phage_name, len(genome))

    if 'tbl' in args.deliverables or len(args.deliverables) == 0:
        # Write features to tbl
        tbl_path = os.path.join(args.output_folder, '%s.tbl' % args.phage_name)
        if test_save(tbl_path):
            write_tbl(df, tbl_path, args.phage_name, len(genome))

    if 'gbf' in args.deliverables and 'sqn' not in args.deliverables:
        raise ValueError("Can not generate gbf without sqn. (Add sqn to desired deliverable list)")

    # Generate sqn file (and gbf if specified)
    if 'sqn' in args.deliverables or len(args.deliverables) == 0:
        from libfaketime import fake_time

        # Check if tbl given too
        if 'tbl' not in args.deliverables and len(args.deliverables) != 0:
            raise ValueError("Can not generate sqn without tbl. (Add tbl to desired deliverable list)")

        if 'fsa' not in args.deliverables and len(args.deliverables) != 0:
            raise ValueError("Can not generate sqn without fsa. (Add fsa to desired deliverable list)")

        # tbl2asn_args = [
        #     # 'tbl2asn',
        #     '/tbltest.sh',
        #     # '-t', SBT_TEMPLATE,   TODO: ACQUIRE TEMPLATE FILE
        #     '-i', fsa_path,
        #     '-j', "[genome=%s][gcode=11]" % args.phage_name,
        #     '-n', args.phage_name,
        #     '-Z', os.path.join(args.output_folder, '%s_tbl2asn_discrepancy_report.txt' % args.phage_name)
        # ]
        #
        # if 'gbf' in args.deliverables or len(args.deliverables) == 0:
        #     tbl2asn_args.append('-V b')

        # import datetime
        # with fake_time('2018-01-01 00:00:01'):
        cmd = 'LD_PRELOAD=/home/daemon/miniconda/envs/mas/lib/python3.8/site-packages/libfaketime/vendor/libfaketime/src/libfaketime.so.1 FAKETIME="2019-01-01 00:00:01" ' \
              f'tbl2asn -i "{fsa_path}" -j [genome="{args.phage_name}"][gcode=11] -n "{args.phage_name}" -Z "{os.path.join(args.output_folder, "%s_tbl2asn_discrepancy_report.txt" % args.phage_name)}" -V b'
        subprocess.run(cmd, check=True, shell=True)

        # Fix locus line
        if 'gbf' in args.deliverables or len(args.deliverables) == 0:
            gbf_path = os.path.join(args.output_folder, '%s.gbf' % args.phage_name)
            first_line = True
            for line in fileinput.input(files=[gbf_path], inplace=True):
                if first_line:
                    line = line.rstrip()
                    after_name_loc = line.find(args.phage_name) + len(args.phage_name)
                    print(line[:after_name_loc] + ' ' + line[after_name_loc:])
                    first_line = False
                else:
                    print(line.rstrip())

#     if 'docx' in args.deliverables or len(args.deliverables) == 0:
#         docx_path = os.path.join(args.output_folder, '%s_Polished.docx' % args.phage_name)
#         if test_save(docx_path):
#             write_docx(df, docx_path, args.phage_name, genome)


if __name__ == '__main__':
    main(parse_args())
