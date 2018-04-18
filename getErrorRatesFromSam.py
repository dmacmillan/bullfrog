__version__ = '1.0'

import argparse
import os
import sys
import re
import logging

output_columns = (
    'qname',
    'qlen',
    'aln_matches',
    'seq_matches',
    'seq_mismatches',
    'insertions',
    'deletions',
    'clips'
)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Given a SAM file, return a list of percent identities ' \
                      'for each alignment based on the cigar string and NM tag. ' \
                      'Warning: this script assumes that given multiple alignments ' \
                      'for the same query that the best alignment is first.',
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
        '-b',
        '--best',
        action = 'store_false',
        help = 'Set this flag to disable selecting only the first alignment per read'
    )

    parser.add_argument(
        '-mcs',
        '--max_clip_start',
        type = int,
        default = '200',
        help = 'Filter any alignments that have a clipped region exceeding this value at the start of the sequence'
    )

    parser.add_argument(
        '-mce',
        '--max_clip_end',
        type = int,
        default = '200',
        help = 'Filter any alignments that have a clipped region exceeding this value at the end of the sequence'
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

    seen = set()

    with open(output_file, 'w') as o:
        o.write(
            ('\t').join(
                [
                    x for x in output_columns
                ]
            ) + '\n'
        )
        with open(args.sam, 'r') as f:
            for i,line in enumerate(f):
                # Skip header lines
                if line[0] == '@':
                    continue
                output_dict = {x:None for x in output_columns}
                cols = line.strip().split('\t')
                qname = cols[0]
                if args.best:
                    if qname in seen:
                        continue
                    seen.add(qname)
                output_dict['qname'] = qname
                flag = int(cols[1])
                cigar = cols[5]
                seq = cols[9]
                qlen = len(seq)
                output_dict['qlen'] = qlen
                # Check if sequence is aligned
                if flag & 4:
                    logging.debug('"{}" is unmapped!'.format(qname))
                    continue
                # Sequence is a secondary alignment
                elif flag & 256:
                    logging.debug('"{}" is a secondary alignment!'.format(qname))
                    continue
                ## 2018-01-10
                # Sequence is a supplementary alignment
                elif flag & 2048:
                    logging.debug('"{}" is a supplementary alignment!'.format(qname))
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
                    logging.warning('"{}" has no NM tag! Continuing'.format(qname))
                    continue
                if cigar == '*':
                    logging.warning('"{}" is mapped, but cigar is "{}". Continuing!'.format(qname, cigar))
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
                clipped_start = clipped_end = 0
                seen_non_clipped = False
                for i,item in enumerate(re.findall(r'\d+[A-Z]', cigar)):
                    op = item[-1]
                    num = int(item[:-1])
                    if (op in ('S', 'H')):
                      if (not seen_non_clipped):
                        clipped_start += num
                      else:
                        clipped_end += num
                    else:
                      seen_non_clipped = i
                    cigar_dict[op] += num
                if (clipped_start > args.max_clip_start):
                  logging.debug('"{}" has clipped start = "{}" > "{}"'.format(qname, clipped_start, args.max_clip_start))
                  continue
                elif (clipped_end > args.max_clip_end):
                  logging.debug('"{}" has clipped end = "{}" > "{}"'.format(qname, clipped_end, args.max_clip_end))
                  continue
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
                output_dict['clips'] = cigar_dict['H'] + cigar_dict['S']
                logging.info('Successfully extracted line {}'.format(i))
                o.write(
                    ('\t').join(
                        [
                            str(output_dict[key]) for key in output_columns
                        ]
                    ) + '\n'
                )
