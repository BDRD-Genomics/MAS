#!/usr/bin/env python
'''
Matt Lueder

This script is for uploading annotations to our internal annotation database. To use this script, create a tsv file
which has 4 columns:
  - Column 1: full amino acid sequence of protein being annotated
  - Column 2: the manual annotation
  - Column 3: notes which we want to make public on genbank
  - Column 4: private notes to be used internally
  - Column 5: flag (green, yellow, red, review name)

To see usage information use the option '-h'.
'''
import sys, os, django

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "LIMS.settings")
django.setup()
import sqlite3
import argparse
import re
import textwrap
import shutil
import math

from amd_database_scripts.paths import *
from amd_database_scripts.flag_category_mapping import *
import create_blast_database


class ComparisonTable:
    '''
    Class for creating a table which compares an annotation that is in the database
    currently to the annotation that will be overwriting it.
    '''
    def __init__(self, seq, currentAnnotation, currentPubNotes, currentPrivNotes, currentFlag,
                 newAnnotation, newPubNotes, newPrivNotes, newFlag):
        '''
        :param seq: Protein sequence
        :param currentAnnotation: The annotation in the database for this sequence
        :param currentPubNotes: The public notes in the database for this sequence
        :param currentPrivNotes: The private notes in the database for this sequence
        :param currentFlag: The flag in the database for this sequence (in text form)
        :param newAnnotation: The annotation being uploaded for this sequence
        :param newPubNotes: The public notes being uploaded for this sequence
        :param newPrivNotes: The private notes being uploaded for this sequence
        :param newFlag: The flag being uploaded for this sequence. (in text form)
        '''
        self.seq = seq
        self.currentAnnotation = currentAnnotation
        self.currentPubNotes = currentPubNotes
        self.currentPrivNotes = currentPrivNotes
        self.currentFlag = currentFlag
        self.newAnnotation = newAnnotation
        self.newPubNotes = newPubNotes
        self.newPrivNotes = newPrivNotes
        self.newFlag = newFlag

        self.headers = ['Field', 'Previous', 'New']

        self.col_width = [len(x) for x in self.headers]
        self.row_height = [1, 1, 1, 1, 1]
        self.rows = [[], [], [], [], []]

    def print_table(self):
        # Determine terminal width
        total_width = shutil.get_terminal_size((80, 80)).columns
        if total_width < 80:
            total_width = 80

        rows = ['Sequence', 'Annotation', 'Public Notes', 'Private Notes', 'Flag']

        # Determine max column width
        col1_size = max([len(x) for x in rows])
        max_col_width = math.floor((total_width - 10 - col1_size) / 2)

        seqF = self.format_sequence(self.seq, max_col_width)
        cAnnotationF = self.format_text(self.currentAnnotation, max_col_width)
        nAnnotationF = self.format_text(self.newAnnotation, max_col_width)
        cPubNotesF = self.format_text(self.currentPubNotes, max_col_width)
        nPubNotesF = self.format_text(self.newPubNotes, max_col_width)
        cPrivNotesF = self.format_text(self.currentPrivNotes, max_col_width)
        nPrivNotesF = self.format_text(self.newPrivNotes, max_col_width)
        cFlagF = self.format_text(self.currentFlag, max_col_width)
        nFlagF = self.format_text(self.newFlag, max_col_width)
        
        self.set_row(0, [[rows[0]], seqF, seqF])
        self.set_row(1, [[rows[1]], cAnnotationF, nAnnotationF])
        self.set_row(2, [[rows[2]], cPubNotesF, nPubNotesF])
        self.set_row(3, [[rows[3]], cPrivNotesF, nPrivNotesF])
        self.set_row(4, [[rows[4]], cFlagF, nFlagF])

        note = 'The protein shown below in the table already has an annotation associated with it.' \
               ' Differences were found in one or more fields.\nIf you would like to keep the ' \
               'previous annotation type "P", or if you would like to keep the new annotation ' \
               'type "N".\nIf you want to select each field individually you can enter a four ' \
               'character sequence where the first character is for the annotation, the second ' \
               'is for the public notes, the third is for the private notes, and the fourth is for ' \
               'the flag. If you choose this option, you can combine the notes fields by selecting "C" ' \
               'instead pf "P" or "N". For example, the character sequence "PNCP" would keep the previous ' \
               'annotation, use the new public notes, combine the private notes (separating them with a ' \
               'semi-colon), and use the previous flag.'
        print('\n'.join(['\n'.join(textwrap.wrap(line, width=total_width, replace_whitespace=False))
                         for line in note.splitlines()]))
        print(self._generate_string())

        # Accept user input
        selected_annotation = None
        selected_pubNote = None
        selected_privNote = None
        selected_flag = None
        while True:
            response = input('[P/N] or [P/N][P/N/C][P/N/C][P/N]:')

            if len(response) != 1 and len(response) != 4:
                print('Your response should be either 1 or 4 characters long.')
                continue

            if response[0] not in ['P', 'N']:
                print('Your response must start with "P" or "N".')
                continue

            if len(response) == 4:
                if response[0] == 'P':
                    selected_annotation = self.currentAnnotation
                elif response[0] == 'N':
                    selected_annotation = self.newAnnotation

                if response[1] == 'P':
                    selected_pubNote = self.currentPubNotes
                elif response[1] == 'N':
                    selected_pubNote = self.newPubNotes
                elif response[1] == 'C':
                    selected_pubNote = '{}; {}'.format(self.currentPubNotes, self.newPubNotes)
                else:
                    print('The second character in your response must be a "P", "N", or "C".')
                    continue

                if response[2] == 'P':
                    selected_privNote = self.currentPrivNotes
                elif response[2] == 'N':
                    selected_privNote = self.newPrivNotes
                elif response[2] == 'C':
                    selected_privNote = '{}; {}'.format(self.currentPrivNotes, self.newPrivNotes)
                else:
                    print('The third character in your response must be a "P", "N", or "C".')
                    continue

                if response[3] == 'P':
                    selected_flag = self.currentFlag
                elif response[3] == 'N':
                    selected_flag = self.newFlag
                else:
                    print('The fourth character in your response must be a "P" or "N".')
                    continue

            else:
                if response == 'P':
                    selected_annotation = self.currentAnnotation
                    selected_pubNote = self.currentPubNotes
                    selected_privNote = self.currentPrivNotes
                    selected_flag = self.currentFlag
                else:
                    selected_annotation = self.newAnnotation
                    selected_pubNote = self.newPubNotes
                    selected_privNote = self.newPrivNotes
                    selected_flag = self.newFlag
            break

        return selected_annotation, selected_pubNote, selected_privNote, selected_flag

    def set_row(self, row_i, values):
        self.rows[row_i] = values

        for col_i in range(len(values)):
            for line in values[col_i]:
                if len(line) > self.col_width[col_i]:
                    self.col_width[col_i] = len(line)

            if len(values[col_i]) > self.row_height[row_i]:
                self.row_height[row_i] = len(values[col_i])

    def _generate_string(self):
        table_string = ''

        def horizontal_rule(line_style):
            line_str = ''
            for i in range(3):
                line_str += '+' + line_style
                for j in range(self.col_width[i]):
                    line_str += line_style
                line_str += line_style
                if i == 2:
                    line_str += '+\n'
            return line_str

        def make_row(values, alignment):
            row_str = ''
            for i in range(len(values)):
                row_str += '| '
                row_str += '{value:{alignment}{width}}'.format(value=values[i],
                                                               width=self.col_width[i],
                                                               alignment=alignment)
                row_str += ' '
                if i == len(values) - 1:
                    row_str += '|\n'
            return row_str

        # Create Header
        table_string += horizontal_rule('=')
        table_string += make_row(self.headers, '^')
        table_string += horizontal_rule('=')

        # Create rows
        for row_i in range(len(self.rows)):
            for line_i in range(self.row_height[row_i]):
                line_vals = []
                for cell in self.rows[row_i]:
                    if len(cell) > line_i:
                        line_vals.append(cell[line_i])
                    else:
                        line_vals.append('')
                table_string += make_row(line_vals, '^')
            table_string += horizontal_rule('-')

        return table_string
                
    @staticmethod
    def format_sequence(seq, max_col_width):
        if len(seq) > max_col_width:
            seq = seq[:max_col_width-6] + ' [...]'

        return [seq]

    @staticmethod
    def format_text(text, max_col_width):
        return textwrap.wrap(text, width=max_col_width)


def parse_args():
    arg_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='This script is for uploading annotations to our internal annotation database.\nTo use this script,'
                    'you can pass the excel from PHATE, or you can create a tsv file which has 4 columns:\n'
                    '  - Column 1: full amino acid sequence of protein being annotated\n'
                    '  - Column 2: the manual annotation\n'
                    '  - Column 3: notes which we want to make public on genbank\n'
                    '  - Column 4: private notes to be used internally\n'
                    '  - Column 5: flag (green, yellow, red, review name)')

    arg_parser.add_argument('annotations', help='PHATE xlsx or a tsv file which has 4 columns - as described above.')

    return arg_parser.parse_args()


def iter_protein(input_file):
    '''
    Iterates through TSV or excel and parses fields
    :param input_file: path to input file
    '''
    if input_file.endswith('.xlsx'):
        import pandas as pd
        df = pd.read_excel(input_file).fillna('')
        for i, row in df.iterrows():
            yield row['Protein Sequence'], row['Final Annotation'], row['Gene Note'], row['Internal Notes'], row['Flag']

    else:
        with open(input_file) as input_fs:
            for line in input_fs:
                sequence, annotation, gbk_notes, private_notes, flag = [x.strip() for x in line.split('\t')]
                yield sequence.upper(), annotation, gbk_notes, private_notes, flag


def is_ascii(string):
    '''
    Determine if a given string is ascii
    :param string: string to test
    :return: true/false
    '''
    return len(string) == len(string.encode())


def main(args):
    '''
    Start script execution
    :param args: Namespace object from argparse
    '''
    # Connect to database
    connection = sqlite3.connect(DB_FILE)
    cursor = connection.cursor()

    # regex patterns for use later in code
    non_base_pattern = re.compile('[^TCGAU]')
    non_aa_pattern = re.compile('[^ACDEFGHIKLMNPQRSTVWY]')

    for sequence, annotation, gbk_notes, private_notes, flag in iter_protein(args.annotations):
        if len(sequence) < 1 or len(annotation) < 1:
            raise RuntimeError('Error: Can not upload sequence or annotation of length zero')

        # Check annotation for illegal characters
        if any([x in annotation for x in ['\n', '|', '\t', '\r', '=', ';']]):
            print("ERROR: The annotation '%s' contains illegal characters. It will not be uploaded." % annotation)
            continue

        # Check if any of the fields contain non-ascii characters
        if any([not is_ascii(x) for x in [sequence, annotation, gbk_notes, private_notes, flag]]):
            print("ERROR: The protein annotated as '%s' contains non-ascii characters in one or more fields."
                  " It will not be uploaded. Please fix this error and re-upload." % annotation)
            continue

        # Check to ensure flag is valid
        try:
            convert_flag_text_to_int(flag)
        except ValueError:
            print("ERROR: The flag '{}' is not valid. The protein with this flag will not be uploaded".format(flag))
            continue

        # If sequence is a protein, ensure correct starting base
        continue_upload = True
        if non_base_pattern.search(sequence) is not None and non_aa_pattern.search(sequence) is None:
            if sequence[0] != 'M':
                print("The sequence %s contains invalid starting amino acid '%s'. Would you like to change this "
                      "amino acid to a 'M'?" % (sequence, sequence[0]))
                while True:
                    response = input('(Y/n) ').lower()
                    if response in ['y', 'yes']:
                        sequence = 'M%s' % sequence[1:]
                        break

                    elif response in ['n', 'no']:
                        print('Would you like to upload this protein anyways? (Not recommended)')
                        while True:
                            response = input('(Y/n) ').lower()
                            if response in ['y', 'yes']:
                                break

                            elif response in ['n', 'no']:
                                continue_upload = False
                                break
                        break

        # If sequence is not a protein, then make sure it is a nucleotide sequence
        elif non_base_pattern.search(sequence) is not None:
            print('The sequence %s contains characters which are not part of the standard unambiguous IUPAC '
                  'DNA, RNA, or Protein alphabet. Would you like to upload anyways?' % sequence)
            while True:
                response = input('(Y/n) ').lower()
                if response in ['y', 'yes']:
                    break

                elif response in ['n', 'no']:
                    continue_upload = False
                    break

        if continue_upload:
            try:
                cursor.execute('INSERT INTO annotations VALUES(?, ?, ?, ?, ?);',
                               (sequence, annotation, gbk_notes, private_notes, convert_flag_text_to_int(flag)))

            except sqlite3.IntegrityError as e:
                # If protein already annotated, we update only if -f arg was supplied
                if "UNIQUE constraint failed: annotations.sequence" in e.__str__():
                    result = cursor.execute('SELECT * FROM annotations WHERE sequence = ?', (sequence,))
                    db_row = result.fetchone()
                    # old_annotation = result.fetchone()[0]

                    # See if there is a difference in the annotation, notes, or flag
                    if annotation != db_row[1] or gbk_notes != db_row[2] or \
                            private_notes != db_row[3] or convert_flag_text_to_int(flag) != db_row[4]:

                        ct = ComparisonTable(sequence, db_row[1], db_row[2], db_row[3],
                                             convert_flag_int_to_text(db_row[4]), annotation, gbk_notes,
                                             private_notes, flag.upper())
                        sel_annotation, sel_pubNotes, sel_privNotes, sel_flag = ct.print_table()
                        cursor.execute(
                            'UPDATE annotations SET annotation = ?, gbk_notes = ?, private_notes = ?, flag = ? WHERE '
                            'sequence = ?;', (sel_annotation, sel_pubNotes, sel_privNotes,
                                              convert_flag_text_to_int(sel_flag), sequence))

                else:
                    print(e)

    connection.commit()
    connection.close()


if __name__ == '__main__':
    main(parse_args())
    # create_blast_db_args = argparse.Namespace(
    #     create_nucl_db=True,
    #     create_protein_db=True,
    #     output_folder=BLAST_DB_LOC
    # )
    # create_blast_database.main(create_blast_db_args)
