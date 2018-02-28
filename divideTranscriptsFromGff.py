import argparse
import os
import sys
import logging
import pysam
from itertools import groupby

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Given a gff file of query sequences mapped to target, ' \
                    'create files for each target containing all mapped ' \
                    'query sequences.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('gff', help='gff file input')
    parser.add_argument('query', help='Fasta file of query transcripts')
    parser.add_argument('target', help='Fasta file of target transcripts')
    parser.add_argument(
        "-l",
        "--log",
        dest="logLevel",
        default='warning',
        choices=[
            'debug',
            'info',
            'warning',
            'error',
            'critical'
        ],
        help='set the logging level'
    )
    parser.add_argument(
        '-o',
        '--outdir',
        default=os.path.join(os.getcwd(), 'divided'),
        help='Path to output to, will create if does not exist'
    )
    
    args = parser.parse_args()
    
    # Create output directory if it does not exist
    if not os.path.isdir(args.outdir):
        try:
            os.makedirs(args.outdir)
        except OSError:
            pass
    
    # Open gff file
    gff_file_object = open(args.gff, 'r')
    
    # Load gff into Pysam object
    gff = pysam.tabix_iterator(
        gff_file_object,
        parser = pysam.asGTF()
    )

    # Load the two fasta files
    query = pysam.FastaFile(args.query)
    target = pysam.FastaFile(args.target)

    # Need to custom split attributes because of missing space in GMAP's output gff
    for key, group in groupby(
        gff,
        key = lambda entry: dict([x.split('=') for x in entry.attributes.split(';')])['Name']
    ):
        transcript = group.next()
        scaffname = transcript.contig
        scaffold_seq_name = os.path.join(args.outdir, '{}.fa'.format(scaffname))
        aligned_seqs_name = os.path.join(args.outdir, 'aligned_to_{}.fa'.format(scaffname))
        if not os.path.isfile(scaffold_seq_name):
            with open(scaffold_seq_name, 'a') as f:
                f.write('>{}\n'.format(scaffname))
                f.write(target.fetch(scaffname) + '\n')
        with open(aligned_seqs_name, 'a') as f:
            f.write('>{}\n'.format(key))
            f.write(query.fetch(key) + '\n')
            
    # Close gff file
    gff_file_object.close()