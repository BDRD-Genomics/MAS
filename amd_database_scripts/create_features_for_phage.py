#!/usr/bin/env python
'''
Matt Lueder

This script creates entries in the feature database table based off results from Glimmer3 and tRNAscan-SE. Annotations
are automatically created from the tRNAscan-SE results, and these are entered into the annotations database (Unless user
specifies otherwise). Features for terminal repeats may also optionally be added. Automated annotations are also created
when terminal repeat features are created.

To see usage information use the option '-h'.
'''
from argparse import Namespace
from amd_database_scripts.genomic_loci_conversions import *


def parse_glimmer_results(path):
    with open(path, 'r') as glimmer_predict:
        for line in glimmer_predict:
            if line.startswith(">"):
                continue
            else:
                lineSplit = [x.strip() for x in line.split()]

                start, stop, strand = glimmer_to_db_standard(
                    start=int(lineSplit[1]),
                    stop=int(lineSplit[2]),
                    strand=lineSplit[3][0]
                )

                yield Namespace(
                    start=start,
                    stop=stop,
                    strand=strand
                )


def parse_trnascan_results(path):
    with open(path, 'r') as trnascan_results:
        ready = False
        for line in trnascan_results:
            if line.startswith('---'):
                ready = True
                continue

            if ready is True:
                lineSplit = [x.strip() for x in line.split()]

                begin = int(lineSplit[2])
                end = int(lineSplit[3])

                # Convert to database standard genomic loci representation
                start, stop, strand = trnascan_to_db_standard(begin, end)

                yield Namespace(
                    start=start,
                    stop=stop,
                    strand=strand,
                    amino=lineSplit[4],
                    codon=lineSplit[5],
                    score=float(lineSplit[8])
                )