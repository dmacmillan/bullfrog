import argparse
import os
import sys
import logging
import re
from customFunctions import *

qi_pattern = re.compile(
    ('').join(
        (
            'QI:',
            '(\d+)\|',
            '([-+]?[0-9]*\.?[0-9]+)\|',
            '([-+]?[0-9]*\.?[0-9]+)\|',
            '([-+]?[0-9]*\.?[0-9]+)\|',
            '([-+]?[0-9]*\.?[0-9]+)\|',
            '([-+]?[0-9]*\.?[0-9]+)\|',
            '(\d+)\|',
            '(\d+)\|',
            '(\d+)'
        )
    )
)

def check_zero(number, mul=1):
    if number == 0:
        return None
    return mul * number

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('fasta', help='fasta file input')
    parser.add_argument("-l", "--log", dest="log_level", default='warning', choices=['debug', 'info', 'warning', 'error', 'critical'], help='set the logging level')
    parser.add_argument(
        '-o',
        '--outdir',
        default=os.path.join(os.getcwd(), 'divided'),
        help='Path to output to, will create if does not exist'
    )
    
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()))
    
    # Create output directory if it does not exist
    if not os.path.isdir(args.outdir):
        try:
            os.makedirs(args.outdir)
        except OSError:
            pass

    for header, seq in iterate_fasta(args.fasta):
        qi = re.search(qi_pattern, header)
        if not qi:
            logging.warning('No QI code for "{}"!'.format(header))
            continue
        utr5 = check_zero(int(qi.group(1)))
        logging.debug('5\'UTR length: {}'.format(utr5))
        utr3 = check_zero(int(qi.group(8)), mul=-1)
        logging.debug('3\'UTR length: {}'.format(utr3))
        print('{}\n{}'.format(header, seq[utr5:utr3]))
        

        
#    with open(args.fasta, 'r') as f:
#        firstline = f.readline().strip()
#        if firstline[0] != '>':
#            logging.critical(
#                'First line of fasta file is not a header!'
#            )
#            sys.exit()
#        seq = ''
#        header = firstline
#        for line in f:
#            if line[0] == '>':
#                qi = re.search(qi_pattern, header)
#                if not qi:
#                    logging.warning('No QI code for "{}"!'.format(header))
#                utr5 = check_zero(int(qi.group(1)))
#                logging.debug('5\'UTR length: {}'.format(utr5))
#                utr3 = check_zero(int(qi.group(8)), mul=-1)
#                logging.debug('3\'UTR length: {}'.format(utr3))
#                if qi:
#                    print '{}\n{}'.format(header, seq[utr5:utr3])
#                seq = ''
#                header = line.strip()
#            else:
#                seq += line.strip()
#        if qi:
#            print '{}\n{}'.format(header, seq[utr5:utr3])
