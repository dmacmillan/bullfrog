import argparse
import os
import sys
import re
import logging
import scipy.stats
import numpy as np

output_columns = (
    'qname',
    'qlen',
    'aln_matches',
    'seq_matches',
    'seq_mismatches',
    'insertions',
    'deletions'
)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Given a SAM file, return a list of percent identities ' \
                      'for each alignment based on the cigar string and NM tag',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('sam', help='SAM file input')
    parser.add_argument(
        '-n',
        '--name',
        help = 'Set the name of the file. Otherwise will be basename of sam ' \
               'appended with ".rates"'
    )
    parser.add_argument(
        '-o',
        '--outdir',
        default = os.getcwd(),
        help = 'Directory to output to, will be created if it does not exist'
    )
    parser.add_argument(
        '-l',
        '--log_level',
        choices=(
            'critical',
            'error',
            'warning',
            'info',
            'debug',
            'notset'
        ),
        default = 'error',
        help='Logging level'
    )
    
    args = parser.parse_args()

    if not args.name:
        args.name = args.sam + '.rates'

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    logging.basicConfig(
        level = getattr(logging, args.log_level.upper())
    )
    
    output_file = os.path.join(args.outdir, args.name)

    with open(output_file, 'w') as o:
        o.write(
            ('\t').join(
                [
                    x for x in output_columns
                ]
            ) + '\n'
        )
        with open(args.sam, 'r') as f:
            for line in f:
                if line[0] == '@':
                    continue
                output_dict = {x:None for x in output_columns}
                cols = line.strip().split('\t')
                qname = cols[0]
                output_dict['qname'] = qname
                flag = int(cols[1])
                cigar = cols[5]
                seq = cols[9]
                qlen = len(seq)
                output_dict['qlen'] = qlen
                # Check if sequence is aligned
                if flag & 4:
                    logging.warning('"{}" is unmapped!'.format(qname))
                    continue
                # Sequence is a secondary alignment
                elif flag & 256:
                    logging.warning('"{}" is a secondary alignment!'.format(qname))
                    continue
                ## 2018-01-10
                # Sequence is a supplementary alignment
                elif flag & 2048:
                    logging.warning('"{}" is a supplementary alignment!'.format(qname))
                    continue
                ##
                # According to my research, the NM tag represents
                # the minimum number of insertions, deletions, and
                # base changes required for the query to become equal
                # To the reference, excluding clipping
                NM = MD = None
                for item in cols[10:]:
                    if item[:2] == 'MD':
                        MD = item[5:]
                    elif item[:2] == 'NM':
                        NM = int(item[5:])
                if not NM:
                    logging.warning('"{}" has no NM tag'.format(qname))
                    continue
                if cigar == '*':
                    logging.warning('"{}" is mapped, but cigar is "{}"'.format(qname, cigar))
                    continue
                cigar_dict = {
                    'M': 0,
                    'I': 0,
                    'D': 0,
                    'N': 0,
                    'S': 0,
                    'H': 0,
                    'P': 0,
                    '=': 0,
                    'X': 0
                }
                for item in re.findall(r'\d+[A-Z]', cigar):
                    op = item[-1]
                    num = int(item[:-1])
                    cigar_dict[op] += num
                output_dict['aln_matches'] = cigar_dict['M']
                # The number of matches is computed as follows:
                output_dict['seq_matches'] = cigar_dict['M'] + \
                           cigar_dict['I'] + \
                           cigar_dict['D'] - \
                           NM
                # Number of mismatches
                output_dict['seq_mismatches'] = cigar_dict['M'] - output_dict['seq_matches']
                output_dict['insertions'] = cigar_dict['I']
                output_dict['deletions'] = cigar_dict['D']
                o.write(
                    ('\t').join(
                        [
                            str(output_dict[key]) for key in output_columns
                        ]
                    ) + '\n'
                )
#            mismatch_rate = nmismatches / float(cigar_dict['M'])
#            insertion_rate = cigar_dict['I'] / float(cigar_dict['M'])
#            deletion_rate = cigar_dict['D'] / float(cigar_dict['M'])

#if not args.read_lengths:
#    for i in range(len(tmat)):
#        print 100 * (tmat[i] + tmis[i] - tnm[i]) / float(tmat[i])
#    mmis = mean(tmis)
#    mins = mean(tins)
#    mdel = mean(tdel)
#    sys.stderr.write('Average mismatch rate: {}\n'.format(mmis))
#    sys.stderr.write('Average insertion rate: {}\n'.format(mins))
#    sys.stderr.write('Average deletion rate: {}\n'.format(mdel))
#    sys.stderr.write('Total average error rate: {}\n'.format(mmis + mins + mdel))
#    
#    #print scipy.stats.describe(tmis)
