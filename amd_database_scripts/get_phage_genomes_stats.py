#!/usr/bin/env python

import sqlite3

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from paths import *
from feature_table import GetFeatureTable


def main():
    '''
    Start script execution
    :param args: Namespace object from argparse
    '''
    # Connect to database
    connection = sqlite3.connect(DB_FILE)
    cursor = connection.cursor()

    # Get all phages
    result = cursor.execute('SELECT * FROM phages')
    phages = result.fetchall()
    if phages is None:
        raise LookupError('ERROR: There are no phages in the database.')

    else:
        for idx, phage in enumerate(phages):
            phages[idx] = [phage[0], Seq(phage[1], IUPAC.unambiguous_dna)]

    polished = []
    unpolished = []
    dtr_phages = []
    for phage in phages:
        df = GetFeatureTable(phage[0], phage[1])
        is_polished = True
        has_dtr = False

        for i, row in df.iterrows():
            # To be polished, all CDS's must have an annotation
            if row['Feature Type'] == 'CDS':
                if not row['Manual Annotation']:
                    is_polished = False

            if row['Feature Type'] == 'repeat_region':
                has_dtr = True

        if is_polished:
            polished.append(phage[0])
        else:
            unpolished.append(phage[0])
        if has_dtr:
            dtr_phages.append(phage[0])

    print('Number of polished in DB = %i' % len(polished))
    print('Number of unpolished in DB = %i' % len(unpolished))
    print('Number of DTR phages in DB = %i\n' % len(dtr_phages))

    print('Polished phages include:')
    for p in polished:
        print(p)
    print()

    print('Unpolished phages include:')
    for p in unpolished:
        print(p)
    print()

    print('Phages containing DTRs include:')
    for p in dtr_phages:
        print(p)
    print()


if __name__ == '__main__':
    main()
