

class HHsearch:
    '''
    Contains full results of hhsearch output.
    '''

    def __init__(self, **kwargs):
        self.query = kwargs['Query']
        self.match_columns = int(kwargs['Match_columns'])
        self.no_of_seqs = kwargs['No_of_seqs']
        self.neff = float(kwargs['Neff'])
        self.searched_hmms = int(kwargs['Searched_HMMs'])
        self.date = kwargs['Date']
        self.command = kwargs['Command']
        self.alignments = []


class HHalignment:
    '''
    An alignment produced by HHsearch (HHpred)
    '''

    def __init__(self, **kwargs):
        self.accession = kwargs['accession']
        self.name = kwargs['name']
        self.prob = float(kwargs['Probab'])
        self.e_value = float(kwargs['E-value'])
        self.score = float(kwargs['Score'])
        self.aligned_cols = int(kwargs['Aligned_cols'])
        self.identities = float(kwargs['Identities'].strip('%')) / 100.0
        self.similarity = float(kwargs['Similarity'])
        self.sum_probs = float(kwargs['Sum_probs'])
        self.template_neff = float(kwargs['Template_Neff'])
        self.tlen = int(kwargs['tlen'])
        self.tstart = int(kwargs['tstart'])
        self.tend = int(kwargs['tend'])
        self.qstart = int(kwargs['qstart'])
        self.qend = int(kwargs['qend'])
        self.text = kwargs['text']

    def __str__(self):
        return self.text


def _parse_alignments(lines):
    '''
    :param lines: List of lines describing the alignment
    :return:
    '''
    start = 0
    arg_dict = {}
    after_gt = False
    after_kv = False
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith('>') and not after_gt:
            arg_dict['accession'] = line[1:line.find(' ')]
            arg_dict['name'] = line[line.find(' ')+1:]
            start = i
            after_gt = True

        elif after_gt:
            for kv in line.split():
                key, value = kv.split('=')
                arg_dict[key] = value
            after_gt = False
            after_kv = True

        elif (line.startswith('No ') or i == len(lines)-1) and after_kv:
            arg_dict['text'] = ''.join(lines[start:i-1])
            yield arg_dict


def _parse_summary_table(lines):
    '''
    Start 'lines' at first row of table (past table header)
    '''
    for line in lines:
        if line == '\n':
            break
        else:
            arg_dict = {}
            ls = line.split()
            arg_dict['tlen'] = int(ls[-1].strip('()'))
            arg_dict['tstart'] = ls[-2].split('-')[0]
            arg_dict['tend'] = ls[-2].split('-')[1]
            arg_dict['qstart'] = ls[-3].split('-')[0]
            arg_dict['qend'] = ls[-3].split('-')[1]
            yield arg_dict


def parse_hhr(file_handle):
    lines = file_handle.readlines()
    file_head = {}
    hh_align_args = []

    for i, line in enumerate(lines):
        # Parse out info in file head
        if line != '\n':
            key, val = line.split(None, 1)
            file_head[key] = val.strip()

        else:
            for d in _parse_summary_table(lines[i+2:]):
                hh_align_args.append(d)
            break

    for i, d in enumerate(_parse_alignments(lines)):
        hh_align_args[i].update(d)

    hhr = HHsearch(**file_head)
    for alignment_arg_set in hh_align_args:
        hhr.alignments.append(HHalignment(**alignment_arg_set))

    return hhr
