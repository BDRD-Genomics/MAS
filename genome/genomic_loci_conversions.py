'''
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


def coordinate_file_to_db_standard(start, stop, strand):
    '''
    Coordinate file rules:
     - 1-based
     - Closed
     - 'strand' represented by '+' or '-'
     - start < stop unless wrap-around
    '''
    if strand == '+':
        return start - 1, stop, strand

    elif strand == '-':
        return stop, start - 1, strand
        # return stop - 1, start - 2, strand

    # elif strand == '+' and start > stop:
    #     return start - 1, stop, strand

    # elif strand == '-' and start > stop:
    #     return stop - 1, start - 2, strand

    else:
        raise ValueError("Invalid coordinates.")


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


def get_protein_sequence(start, stop, strand, genome, table=11):
    '''
    Extract protein sequence from genome given coordinates in db standard format.
    'genome' is a Bio.Seq object
    '''
    #  Not using cds=true because it causes issues with ambiguous seqs
    # return get_dna_sequence(start, stop, strand, genome).translate(table=11, cds=True)
    prot = get_dna_sequence(start, stop, strand, genome).translate(table=table, to_stop=True)
    return 'M' + prot[1:]


def get_rna_sequence(start, stop, strand, genome):
    return get_dna_sequence(start, stop, strand, genome).transcribe()